import os
import json
import requests
from tqdm import tqdm
import xml.etree.ElementTree as ET
from tenacity import retry, stop_after_attempt, wait_fixed, retry_if_exception_type


def check_item(item):
    if item is None:
        return "NaN"
    else:
        return item.text


# 使用 tenacity 的装饰器来定义重试策略：
@retry(
    stop=stop_after_attempt(6),  # 重试 6 次
    wait=wait_fixed(3),  # 每次重试之间等待 3 秒
    retry=retry_if_exception_type(
        requests.RequestException
    ),  # 当发生 requests.RequestException 类型的异常时就重试
)
def call_efetch(gene_id):
    """调用 NCBI Web APIs 的 EFetch 方法, 根据输入的 Entrez Gene ID 获取基因相关信息。
    官方文档: https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_

    输入：
      gene_id: Entrez Gene ID

    输出：
      response: 根据输入的 Gene ID 获得基因相关信息的响应

    """

    # EFetch 方法的 Base URL
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    my_api_key = "your_api_key********************************"

    """参数说明

    db: 要进行搜索的数据库 gene, 值必须是有效的 Entrez 数据库名称。https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    id: UID 列表。可以提供单个 UID 或以逗号分隔的 UIDs 列表。
    retmode: 检索类型, 确定返回输出的格式。这里使用 XML 格式返回结果。

    示例: Fetch full XML record for Gene ID 2
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=2&retmode=xml
    """
    params = {"db": "gene", "id": gene_id, "retmode": "xml", "api_key": my_api_key}
    response = requests.get(base_url, params=params)
    # print(response.url)
    response.raise_for_status()  # 检查响应对象中是否发生错误

    return response


os.makedirs("./genes_xml", exist_ok=True)
os.makedirs("./logs", exist_ok=True)
# 基因符号转换到 Entrez Gene ID 的字典
with open("./data/gene_to_id.json", "r") as f:
    gene_id_dict = json.load(f)

# 记录程序运行的日志和初始化存储结果的字典：
log_file = open("./logs/genes_info_20240209.txt", "a+")
log_file.write(f"基因数量: {len(gene_id_dict.keys())}" + "\n\n")
gene_info_dict = {}
nums_limit = 100000
num = 1

for gene_name, gene_id in tqdm(gene_id_dict.items()):
    if num == nums_limit:
        break

    try:
        response = call_efetch(gene_id)
        gene_xml = response.text
        with open(f"./genes_xml/{gene_name}__{gene_id}.xml", "w") as f:
            f.write(gene_xml)

        root = ET.fromstring(gene_xml)
        # 官方符号通常是唯一的：
        official_symbol = root.find(".//Entrezgene_gene/Gene-ref/Gene-ref_locus").text
        # 充分考虑到基因别名：
        gene_aliases = root.findall(
            ".//Entrezgene_gene/Gene-ref/Gene-ref_syn/Gene-ref_syn_E"
        )
        if gene_aliases == []:
            gene_aliases = []
        else:
            gene_aliases = [g.text for g in gene_aliases]
        # 基因类型：
        gene_type = root.find(".//Entrezgene_type")
        if gene_type is None:
            gene_type = "NaN"
        else:
            gene_type = gene_type.attrib["value"]  # 取属性的 value 值
        # 基因的正式描述：
        description = check_item(root.find(".//Entrezgene_gene/Gene-ref/Gene-ref_desc"))
        # 基因的 NCBI Summary 信息：
        gene_summary = check_item(root.find(".//Entrezgene_summary"))

        gene_info_dict[gene_id] = {
            "official_symbol": official_symbol,
            "description": description,
            "gene_type": gene_type,
            "summary_info": gene_summary,
            "gene_aliases": gene_aliases,
        }

        log_file.write(
            f"对于基因 {gene_name}, Gene ID: {gene_id}, 成功获取到关键信息!" + "\n\n"
        )

    except Exception as e:
        log_file.write(
            f"在尝试 6 次后，基因 {gene_name} | Gene ID 为 {gene_id}, 查询相关信息失败。错误：{e}"
            + "\n\n"
        )
        print(response.text)
        with open("failed_genes_id.txt", "a+") as f:
            f.write(f"{gene_name}__{gene_id}" + "\n")

    num += 1

# 关闭文件对象：
log_file.close()
# 保存嵌套字典数据：
json_str = json.dumps(gene_info_dict, ensure_ascii=False, indent=2)
with open("./data/genes_info.json", "w") as f:
    f.write(json_str)
