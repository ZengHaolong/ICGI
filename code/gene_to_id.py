import json
import requests
from tqdm import tqdm
import xml.etree.ElementTree as ET
from tenacity import retry, stop_after_attempt, wait_fixed, retry_if_exception_type


# 使用 tenacity 的装饰器来定义重试策略：
@retry(
    stop=stop_after_attempt(6),  # 重试 6 次
    wait=wait_fixed(3),  # 每次重试之间等待 3 秒
    retry=retry_if_exception_type(
        requests.RequestException
    ),  # 当发生 requests.RequestException 类型的异常时就重试
)
def call_esearch(gene_name, max_uids, sort_uids):
    """调用 NCBI Web APIs 的 ESearch 方法, 根据输入的基因符号查询到准确的 Entrez Gene ID。
    官方文档: https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_

    输入：
      gene_name: 基因符号
      max_uids: 从检索集中检索到的 UIDs 的最大数量
      sort_uids: 在 ESearch 输出中用于排序 UID 的方法

    输出：
      response: 根据输入的基因符号获得 Gene ID 的响应

    """

    # ESearch 方法的 Base URL
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    my_api_key = "your_api_key********************************"

    """参数说明

    db: 要进行搜索的数据库 gene。值必须是有效的 Entrez 数据库名称。https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    term: 输入 Entrez 文本查询。所有特殊字符必须进行 URL 编码。
    retmode: 检索类型。决定了返回输出的格式。这里使用 JSON 格式返回结果。
    retmax: 从检索集中检索到的 UIDs 的最大数量 (默认为 20)。这里设置为 25。
    sort: 指定了在 ESearch 输出中对 UIDs 进行排序的方法。这里使用 relevance。

    限制 term 的搜索字段。如果使用, 整个搜索词将限制在指定的 Entrez 字段中。
    "term": f"{gene_name}[GENE] AND Homo sapiens[ORGN]",
    "term": f"{gene_name}[Gene Name] AND Homo sapiens[Organism]",
    """
    params = {
        "db": "gene",
        "retmode": "json",
        "retmax": max_uids,
        "term": f"{gene_name}[GENE] AND Homo sapiens[ORGN]",
        "sort": sort_uids,
        "api_key": my_api_key,
    }
    response = requests.get(base_url, params=params)
    # print(response.url)
    response.raise_for_status()  # 检查响应对象中是否发生错误

    return response


@retry(
    stop=stop_after_attempt(6),
    wait=wait_fixed(3),
    retry=retry_if_exception_type(requests.RequestException),
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


with open("./data/query_genes.txt", "r") as f:
    all_genes = f.read().split("\n")

# 记录程序运行的日志：
log_file = open("./logs/log_20240208.txt", "a+")
log_file.write(f"基因集的数量: {len(all_genes)}" + "\n\n")

# 基因符号转换到 Gene ID 的字典
gene_id_dict = {}

for gene in tqdm(all_genes):
    try:
        response1 = call_esearch(gene, max_uids=25, sort_uids="relevance")
        # print(response1.json())
        ids_list = response1.json()["esearchresult"]["idlist"]
        candidate_gene_id = []

        if len(ids_list) == 0:
            log_file.write(
                f"对于基因 {gene}, 未在 NCBI 中查询到准确的 Gene ID! " + "\n\n"
            )
            with open("genes_None.txt", "a+") as f:
                f.write(gene + "\n")
            continue  # continue 语句用于跳过当前循环迭代中的剩余部分, 并立即开始下一次迭代。

        elif len(ids_list) == 1:
            gene_id = ids_list[0]

            response2 = call_efetch(gene_id)
            root = ET.fromstring(response2.text)
            discontinued_flag = root.find(
                ".//Entrezgene_track-info/Gene-track/Gene-track_discontinue-date"
            )
            # 在使用 xml.etree.ElementTree 的 find() 方法时没有找到指定的元素，该方法会返回 None。
            if discontinued_flag is None:  # 正常
                log_file.write(
                    f"对于基因 {gene}, 在 NCBI 中查询到唯一的 Gene ID: {ids_list[0]}"
                    + "\n\n"
                )
            else:  # 关注到可能已经被 NCBI Discontinued 的 Gene ID
                log_file.write(
                    f"查询到基因 {gene} 的唯一 Gene ID {check_gene_id}, 但被 NCBI Discontinue 了!"
                    + "\n\n"
                )
                with open("genes_NCBI_Discontinued_1.txt", "a+") as f:
                    f.write(f"{gene}__{gene_id}" + "\n")
                continue  # continue 语句用于跳过当前循环迭代中的剩余部分, 并立即开始下一次迭代。

        elif len(ids_list) > 1:
            log_file.write(
                f"对于基因 {gene}, 在 NCBI 中查询到 {len(ids_list)} 个 Gene ID: {ids_list}"
                + "\n\n"
            )
            gene_id = ids_list[0]

            for check_gene_id in ids_list:
                response2 = call_efetch(check_gene_id)
                root = ET.fromstring(response2.text)

                discontinued_flag = root.find(
                    ".//Entrezgene_track-info/Gene-track/Gene-track_discontinue-date"
                )
                if discontinued_flag is None:
                    log_file.write(
                        f"检查对于基因 {gene} 的 Gene ID {check_gene_id}..." + "\n\n"
                    )
                else:  # 不考虑已经被 NCBI Discontinued 的 Gene ID
                    log_file.write(
                        f"对于基因 {gene} 的 Gene ID {check_gene_id}, 已经被 NCBI Discontinue 了!"
                        + "\n\n"
                    )
                    with open("genes_NCBI_Discontinued_2.txt", "a+") as f:
                        f.write(f"{gene}__{check_gene_id}" + "\n")
                    continue

                # 在使用 xml.etree.ElementTree 的 findall() 方法时，如果没有找到指定的元素，该方法会返回一个空列表 []。
                # 官方符号通常是唯一的
                official_symbol = root.findall(
                    ".//Entrezgene_gene/Gene-ref/Gene-ref_locus"
                )
                # 充分考虑到基因别名（许多基因具有多个别名）
                gene_aliases = root.findall(
                    ".//Entrezgene_gene/Gene-ref/Gene-ref_syn/Gene-ref_syn_E"
                )

                if gene_aliases == []:
                    gene_aliases = {}
                else:
                    gene_aliases = {g.text for g in gene_aliases}

                if official_symbol == []:
                    official_symbol = []
                else:
                    official_symbol = [o.text for o in official_symbol]

                if not official_symbol:
                    log_file.write(
                        f"对于基因 {gene} 的 Gene ID {check_gene_id}, 未查询到基因的 Official Symbol!"
                        + "\n\n"
                    )
                elif len(official_symbol) > 1:
                    log_file.write(
                        f"对于基因 {gene} 的 Gene ID {check_gene_id}, Official Symbol 有多个!"
                        + "\n\n"
                    )
                    official_symbol = official_symbol[0]
                else:
                    official_symbol = official_symbol[0]

                if gene == official_symbol:  # 基因符号与 Official Symbol 匹配最重要
                    candidate_gene_id = []  # candidate_gene_id 就置为空列表
                    gene_id = check_gene_id
                    break  # 中断内层循环 -> for check_gene_id in ids_list:
                else:
                    if gene in gene_aliases:
                        candidate_gene_id.append(check_gene_id)

        if (
            len(candidate_gene_id) > 0
        ):  # 若未找到与 Official Symbol 匹配的基因名，就考虑别名匹配，选择 relevance 排最前面的。
            gene_id = candidate_gene_id[0]

        # 基因符号: Gene ID
        gene_id_dict[gene] = gene_id
        log_file.write(
            f"对于基因 {gene} 的 Entrez Gene ID {gene_id}, 转换为 Gene ID 成功!"
            + "\n\n"
        )

    except Exception as e:
        log_file.write(
            f"在尝试 6 次后，基因 {gene} 转换为 Gene ID 的转化失败。错误：{e}" + "\n\n"
        )
        with open("genes_None.txt", "a+") as f:
            f.write(gene + "\n")

# 关闭文件对象：
log_file.close()

json_str = json.dumps(gene_id_dict, ensure_ascii=False, indent=2)
with open("./data/gene_to_id.json", "w") as f:
    f.write(json_str)
