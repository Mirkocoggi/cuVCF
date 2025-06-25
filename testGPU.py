# source rapids_env/bin/activate
import cupy as cp
import GPUParser as vcf
import cudf as pd
import time
import numpy as np
import math
import gc

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


res = vcf.vcf_parsed()

print("Start parsing")

start_run = time.perf_counter()
res.run("data/felis_catus.vcf", 16)
end_run = time.perf_counter()
run_time = end_run - start_run
print(f"Tempo di esecuzione di run: {run_time:.4f} secondi")

start_df = time.perf_counter()
data1 = vcf.get_var_columns_data(res.var_columns)
data2 = vcf.get_alt_columns_data(res.alt_columns)
data3 = vcf.get_sample_columns_data(res.samp_columns)
data4 = vcf.get_alt_format_data(res.alt_sample)
end_df = time.perf_counter()
df_time = end_df - start_df
print(f"Tempo di creazione dei DataFrame c++: {df_time:.4f} secondi")

start_dfpd = time.perf_counter()
df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
end_dfpd = time.perf_counter()
dfpd_time = end_dfpd - start_dfpd
print(f"Tempo di creazione dei DataFrame pandas: {dfpd_time:.4f} secondi")

start_csv = time.perf_counter()
npartitions = 10  # puoi regolare questo parametro in base alle tue necessità
save_cudf_to_csv_in_chunks(df1, "df1.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df2, "df2.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df3, "df3.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df4, "df4.csv", npartitions, index=False)
end_csv = time.perf_counter()
csv_time = end_csv - start_csv
print(f"Tempo di salvataggio in CSV: {csv_time:.4f} secondi")

print("CSV creati")
del df1
del df2
del df3
del df4
gc.collect()
print("Dataframe eliminati")

time_elapsed = 0.0
start_run = time.perf_counter()
df1 = pd.read_csv("df1.csv", delimiter=",")
df2 = pd.read_csv("df2.csv", delimiter=",")
df3 = pd.read_csv("df3.csv", delimiter=",")
df4 = pd.read_csv("df4.csv", delimiter=",")
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)

print(f"Tempo caricamento da csv: {time_elapsed:.4f} secondi")