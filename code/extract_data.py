import json
import pandas as pd
from tqdm import tqdm
from collections import Counter


def is_zero(num, tolerance=1e-6):
    if isinstance(num, int):  # Integer check
        return num == 0
    # 对于浮点数，判断其绝对值是否小于给定的容差（默认为 1e-6）。
    elif isinstance(num, float):  # Float check
        return abs(num) < tolerance
    else:
        raise ValueError("Unsupported numeric type")


def load_metadata(dataset):
    """加载指定数据集的元数据文件（JSON 格式）"""
    with open(f"./{dataset}/metadata.cart.2023-12-25.json", "r") as f:
        return json.load(f)


def load_gene_ids():
    """加载基因 ID 数据，并生成一个从基因 ID 到基因名称的映射字典。"""
    with open("gene_ensg_id_dict.json", "r") as f:
        gene_id_data = json.load(f)
    with open("gene_to_id.json", "r") as f:
        tmp_gene_id = json.load(f)
    selected_genes = {i for i in tmp_gene_id.keys()}
    gene_id_data = [(k, v) for k, v in gene_id_data.items() if k in selected_genes]
    id_to_gene = {v: k for k, v in gene_id_data}
    gene_ids = [gene_and_id[1] for gene_and_id in gene_id_data]

    return gene_ids, id_to_gene


# 处理 TCGA 项目的 RNA-Seq 数据集，提取基因表达矩阵和样本标签，并将处理后的数据保存为 CSV 文件。
for dataset in ["LUAD", "LUSC", "BLCA", "BRCA", "KIRC", "LIHC"]:
    print("-" * 100)
    print(f"Dataset: {dataset}")

    # 加载元数据并提取样本的 TCGA ID 和 Case ID，统计这些 ID 的数量。
    cart = load_metadata(dataset)
    tmp1 = {item["associated_entities"][0]["entity_submitter_id"] for item in cart}
    tmp2 = {item["associated_entities"][0]["case_id"] for item in cart}
    print(f"TCGA ID Nums: {len(tmp1)} \t Case ID Nums: {len(tmp2)}")
    tmp3 = [
        "-".join(item["associated_entities"][0]["entity_submitter_id"].split("-")[:4])
        for item in cart
    ]
    print(Counter(tmp3).most_common(10))

    samples, data = [], []
    samples_set = set()
    for item in tqdm(cart):  # 遍历字典列表
        file_id = item["file_id"]  # File ID
        file_name = item["file_name"]  # TSV Filename
        tcga_id = item["associated_entities"][0]["entity_submitter_id"]  # TCGA ID
        case_id = item["associated_entities"][0]["case_id"]  # Case ID

        sample_vial = tcga_id.split("-")[3]  # Sample-Vial
        # 01A: Primary Solid Tumor      11A: Solid Tissue Normal
        if sample_vial not in {"01A", "11A"}:
            continue
        sample_id = "-".join(tcga_id.split("-")[:4])
        if sample_id in samples_set:
            continue

        sample_path = f"./{dataset}/samples_info/{file_id}/{file_name}"
        df = pd.read_csv(sample_path, sep="\t", skiprows=1)
        df = df.iloc[4:, :]
        df = df[~df["gene_id"].str.contains("_PAR_Y")]  # 删除重复基因
        df.reset_index(drop=True, inplace=True)

        sample_label = 1 if sample_vial == "01A" else 0
        # ENSG Gene ID 与对应的基因表达量 TPM
        tmp_dict = {
            ensg_id: gene_expression
            for ensg_id, gene_expression in zip(
                df["gene_id"], df["tpm_unstranded"]
            )
        }
        sample_info = {
            k.split(".")[0]: tmp_dict[k] for k in sorted(tmp_dict.keys())
        }
        sample_info["label"] = sample_label

        samples.append(tcga_id)
        samples_set.add(sample_id)
        data.append(sample_info)

    # 通过字典列表来创建 DataFrame
    data_df = pd.DataFrame(data, index=samples)
    print(data_df.iloc[:, :-1].shape)

    gene_ids, id_to_gene = load_gene_ids()
    # 只保留能在 NCBI Gene 数据库中通过基因名称准确查询到 Entrez Gene ID 的基因。
    # 根据基因集提取基因表达矩阵：
    data_df1 = data_df[gene_ids + ["label"]].copy(deep=True)  # 深拷贝提取基因表达数据
    data_df1.columns = [id_to_gene[ensg] for ensg in gene_ids] + ["label"]
    print(data_df1.iloc[:, :-1].shape)

    # 过滤低表达基因：排除那些表达量为 0 的样本数占总样本数比例超过 50% 的基因。
    genes_selected = []
    N, _ = data_df1.shape
    for col in data_df1.iloc[:, :-1].columns:
        # percentage_0 = (data_df1[col] == 0).sum() / N
        # 计算该列中值为 0 的元素的个数占总行数的比例
        percentage_0 = sum([1 if is_zero(value) else 0 for value in data_df1[col]]) / N
        # 如果这个比例大于 0.50，那么就跳过这一列，否则就将这一列的名称添加到 genes_selected 列表中。
        # if percentage_0 > 0.50 or data_df1[col].mean() < 1:
        if percentage_0 > 0.50:
            continue
        else:
            genes_selected.append(col)
    genes_selected.append("label")

    # genes_selected = [col for col in data_df1.columns if data_df1[col].mean() >= 1 and col != "label"]
    # genes_selected.append("label")

    # 基因表达矩阵与样本标签（癌症/正常）共同用于后续的分析和建模。
    data_df1 = data_df1[genes_selected]
    y = data_df1["label"]
    print("After processing:", data_df1.iloc[:, :-1].shape)
    print(y.value_counts())
    # 将处理后的基因表达矩阵和样本标签保存为 CSV 文件：
    data_df1.to_csv(f"./dataset/{dataset}_TPM.csv", encoding="utf-8")
