import pandas as pd
import matplotlib.pyplot as plt


def get_number(rank_code):
    return int(rank_code[1:]) if len(rank_code) > 1 else 0


def get_best_ranks(df):
    df = df.query('rank_code.str.contains("S")', engine='python')
    df.reset_index(drop=True, inplace=True)

    idx = 0
    df_filtered = pd.DataFrame()

    while idx < len(df.index) - 1:
        row = df.iloc[[idx]]
        next_row = df.iloc[[idx + 1]]

        if get_number(row["rank_code"].item()) >= get_number(next_row["rank_code"].item()):
            df_filtered = pd.concat([df_filtered, row], ignore_index=True)
        idx += 1
    df_filtered = pd.concat([df_filtered, df.iloc[[idx]]], ignore_index=True)

    return df_filtered


def recalculate_percentages(df):
    total = df["fragments_assigned"].sum()
    df["percentage"] = df["fragments_assigned"] / total
    return df


def get_species_df(file_path):
    data = pd.read_csv(file_path, sep="\t", header=None)

    data.columns = ['percentage', 'fragments_covered', 'fragments_assigned', "rank_code", "taxon_id", "taxon_name"]

    # only species
    data = get_best_ranks(data)
    data = recalculate_percentages(data)
    data = data.sort_values('percentage', ascending=False)

    return data


def draw_plot(df, plot_file_path):
    plt.figure(figsize=(12, 8))
    plt.bar(df['taxon_name'].astype(str), df["percentage"], color="skyblue", log=True)
    plt.xlabel("Taxon ID")
    plt.ylabel("Percent")
    plt.title(f"Percentage Distribution of Taxon IDs")
    plt.xticks(rotation=90)
    plt.ylim(0, max(df["percentage"]) + 0.01)
    plt.tight_layout()

    plt.savefig(plot_file_path)


def analyze_kraken_output(file_path, plot_path, filtered_taxon_ids_path, min_percent):
    data = get_species_df(file_path)
    df_filtered = data.query(f"percentage > {min_percent}")

    draw_plot(df_filtered, plot_path)

    with open(filtered_taxon_ids_path, 'w') as output:
        for taxon_id in df_filtered["taxon_id"]:
            output.write(f"{taxon_id}\n")


def parse(input_file, output_folder, taxon_ids, i):
    import os
    import re

    # Создание выходной папки, если она не существует
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Шаблон для поиска таксономического идентификатора
    taxon_pattern = re.compile(r"kraken:taxid\|(\d+)")
    with open(taxon_ids, 'r') as file:
        taxon_list = [line.strip() for line in file]
    # Словарь для хранения данных по таксонам
    taxon_data = {}
    with open(input_file, 'r') as file:
        while True:
            # Чтение блока из 4 строк
            id_line = file.readline().strip()
            if not id_line:
                break  # Конец файла
            seq_line = file.readline().strip()
            plus_line = file.readline().strip()
            qual_line = file.readline().strip()

            # Извлечение таксономического идентификатора
            match = taxon_pattern.search(id_line)
            if match:
                taxon_id = match.group(1)
                if taxon_id in taxon_list:

                    if taxon_id not in taxon_data:
                        taxon_data[taxon_id] = []
                    # Сохранение блока данных
                    taxon_data[taxon_id].append((id_line, seq_line, plus_line, qual_line))

    # Запись данных по таксонам в отдельные файлы
    for taxon_id, records in taxon_data.items():
        output_dir = f"{output_folder}/taxon_{taxon_id}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_file = os.path.join(output_dir, f"taxon_{taxon_id}_{i}.fq")
        with open(output_file, 'w') as out_f:
            for record in records:
                out_f.write("\n".join(record) + "\n")


def parse_fq_by_taxon(input_1, input_2, taxon_id, out_dir):
    parse(input_1, out_dir, taxon_id, 1)
    parse(input_2, out_dir, taxon_id, 2)