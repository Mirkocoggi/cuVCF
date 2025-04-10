# source rapids_env/bin/activate

import GPUParser as vcf
import cudf
import time
import numpy as np
import math

def save_cudf_to_csv_in_chunks(df, filename, npartitions=10, index=False):
    """
    Salva il DataFrame cuDF in un unico file CSV, scrivendo i dati a chunk.
    
    Parameters:
        df (cudf.DataFrame): DataFrame da salvare.
        filename (str): Nome del file CSV finale.
        npartitions (int): Numero di chunk in cui suddividere il DataFrame.
        index (bool): Se salvare o meno l'indice nel CSV.
    """
    n_rows = len(df)
    chunk_size = math.ceil(n_rows / npartitions)
    
    # Ottieni l'intestazione per scriverla solo una volta
    header = ",".join(df.columns.astype(str)) + "\n"
    
    with open(filename, "w") as f:
        f.write(header)
        for i in range(npartitions):
            start = i * chunk_size
            end = min(start + chunk_size, n_rows)
            chunk = df.iloc[start:end]
            # Converti il chunk in CSV senza header
            csv_str = chunk.to_csv(index=index, header=False)
            f.write(csv_str)

# Parsing e misurazione dei tempi
res = vcf.vcf_parsed()
print("Start parsing")
start_run = time.perf_counter()
res.run("data/felis_catus.vcf", 16)
end_run = time.perf_counter()
print(f"Tempo di esecuzione di run: {end_run - start_run:.4f} secondi")

start_df = time.perf_counter()
data1 = vcf.get_var_columns_data(res.var_columns)
data2 = vcf.get_alt_columns_data(res.alt_columns)
data3 = vcf.get_sample_columns_data(res.samp_columns)
data4 = vcf.get_alt_format_data(res.alt_sample)
end_df = time.perf_counter()
print(f"Tempo di creazione dei dizionari in C++: {end_df - start_df:.4f} secondi")

# Converti "qual" da uint16 a float32, se necessario
if "qual" in data1:
    arr = data1["qual"].view(np.float16)
    data1["qual"] = arr.astype(np.float32)

start_dfcudf = time.perf_counter()
df1 = cudf.DataFrame(data1)
df2 = cudf.DataFrame(data2)
df3 = cudf.DataFrame(data3)
df4 = cudf.DataFrame(data4)
end_dfcudf = time.perf_counter()
print(f"Tempo di creazione dei DataFrame cuDF: {end_dfcudf - start_dfcudf:.4f} secondi")

# Salva i DataFrame in un unico file CSV per ciascuno, suddividendo internamente in chunk
start_csv = time.perf_counter()

npartitions = 10  # puoi regolare questo parametro in base alle tue necessità
save_cudf_to_csv_in_chunks(df1, "df1.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df2, "df2.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df3, "df3.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df4, "df4.csv", npartitions, index=False)

end_csv = time.perf_counter()
print(f"Tempo di salvataggio in CSV: {end_csv - start_csv:.4f} secondi")
