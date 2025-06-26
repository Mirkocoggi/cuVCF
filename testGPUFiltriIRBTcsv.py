# source rapids_env/bin/activate
import cupy as cp
import GPUParser as vcf
import cudf
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

time_elapsed = 0

# Parsing e misurazione dei tempi
res = vcf.vcf_parsed()
print("Start parsing")
start_run = time.perf_counter()
res.run("data/IRBT.vcf", 16)
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)

print(f"Parsing time: {time_elapsed:.4f} secondi")

start_run = time.perf_counter()
data1 = vcf.get_var_columns_data(res.var_columns) 
data2 = vcf.get_alt_columns_data(res.alt_columns)
data3 = vcf.get_sample_columns_data(res.samp_columns)
data4 = vcf.get_alt_format_data(res.alt_sample)


# Converti "qual" da uint16 a float32, se necessario
if "qual" in data1:
    arr = data1["qual"].view(np.float16)
    data1["qual"] = arr.astype(np.float32)

# patch perchè non va il binding su var_id
n = len(data3["var_id"])
group_size = 8
data3["var_id"] = np.repeat(np.arange((n + group_size - 1) // group_size), group_size)[:n]

df1 = cudf.DataFrame(data1)
df2 = cudf.DataFrame(data2)
df3 = cudf.DataFrame(data3)
df4 = cudf.DataFrame(data4)

end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Dataframe creati in: {time_elapsed:.4f} secondi")

start_run = time.perf_counter()
npartitions = 10  # puoi regolare questo parametro in base alle tue necessità
save_cudf_to_csv_in_chunks(df1, "df1.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df2, "df2.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df3, "df3.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df4, "df4.csv", npartitions, index=False)
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)

print(f"CSV creati in: {time_elapsed:.4f} secondi")

del df1
del df2
del df3
del df4

gc.collect()
print("Dataframe eliminati")

filter_map = vcf.GTMapGlobal
ref_value = filter_map["1|1"]

time_elapsed = 0.0

start_run = time.perf_counter()
df1 = cudf.read_csv("df1.csv", delimiter=",")
df2 = cudf.read_csv("df2.csv", delimiter=",")
df3 = cudf.read_csv("df3.csv", delimiter=",")
df4 = cudf.read_csv("df4.csv", delimiter=",")
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)

print(f"Tempo caricamento da csv: {time_elapsed:.4f} secondi")

df2["AF"] = df2["AF"].astype(np.uint16)

tot_time = 0
#Filter AC>8
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo AC>8: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed


del df22
gc.collect()

#Filter AC > 8 & AF > 0.5
af_cupy = cp.asarray(df2["AF"].values)
af_float16 = af_cupy.view(cp.float16)
mask = af_float16 > cp.float16(0.5)
start_run = time.perf_counter()
df22 = df2[mask.get()]
df22 = df22[df22["AC"] > 8]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo AC > 8 and AF > 0.5: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed


del df22
gc.collect()

#Filter GT = 1|1
start_run = time.perf_counter()
df33 = df3[df3["GT0"] == ref_value]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo GT = 1|1: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
gc.collect()

#Filter AD[0] > 3
start_run = time.perf_counter()
df33 = df3[df3["AD0"] > 3]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AD[0] > 3: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
gc.collect()

#Filter AD[0] > 3 & AD[1] < 10
start_run = time.perf_counter()
df33 = df3[(df3["AD0"] > 3) & (df3["AD1"] < 10)]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AD[0] > 3 & AD[1] < 10: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
gc.collect()

#Filter GT = 1|1 & AC>8
start_run = time.perf_counter()
df33 = df3[df3["GT0"] == ref_value]
df22 = df2[df2["AC"] > 8]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  GT = 1|1 & AC>8: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter GT = 1|1 & AC>8 & AF > 0.5
af_cupy = cp.asarray(df2["AF"].values)
af_float16 = af_cupy.view(cp.float16)
mask = af_float16 > cp.float16(0.5)
start_run = time.perf_counter()
df33 = df3[df3["GT0"] == ref_value]
df22 = df2[mask.get()]
df22 = df22[df22["AC"] > 8]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  GT = 1|1 & AC>8 & AF > 0.5: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AD[0] > 3 & AC>8
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[df3["AD0"] > 3]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AD[0] > 3 & AC>8: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AC > 8 & AD[0] > 3 & AD[1] < 10
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[(df3["AD0"] > 3) & (df3["AD1"] < 10)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & AD[0] > 3 & AD[1] < 10: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AC > 8 & AF > 0.5 & AD[0] > 3
af_cupy = cp.asarray(df2["AF"].values)
af_float16 = af_cupy.view(cp.float16)
mask = af_float16 > cp.float16(0.5)
start_run = time.perf_counter()
df22 = df2[mask.get()]
df22 = df22[df22["AC"] > 8]
df33 = df3[df3["AD0"] > 3]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & AF > 0.5 & AD[0] > 3: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AC > 8 & AF > 0.5 & AD[0] > 3 & AD[1] < 10
af_cupy = cp.asarray(df2["AF"].values)
af_float16 = af_cupy.view(cp.float16)
mask = af_float16 > cp.float16(0.5)
start_run = time.perf_counter()
df22 = df2[mask.get()]
df22 = df22[df22["AC"] > 8]
df33 = df3[(df3["AD0"] > 3) & (df3["AD1"] < 10)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & AF > 0.5 & AD[0] > 3 & AD[1] < 10: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AD[0] > 3 & GT = 1|1
start_run = time.perf_counter()
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3)]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AD[0] > 3 & GT = 1|1: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
gc.collect()

#Filter AD[0] > 3 & AD[1] < 10 & GT = 1|1
start_run = time.perf_counter()
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3) & (df3["AD1"] < 10)]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AD[0] > 3 & AD[1] < 10 & GT = 1|1: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
gc.collect()

#Filter AC > 8 & GT = 1|1 & AD[0] > 3
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & GT = 1|1 & AD[0] > 3: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AC > 8 & AD[0] > 3 & AD[1] < 10 & GT = 1|1
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3) & (df3["AD1"] < 10)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & AD[0] > 3 & AD[1] < 10 & GT = 1|1: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AC > 8 & AF > 0.5 & GT = 1|1 & AD > 3
af_cupy = cp.asarray(df2["AF"].values)
af_float16 = af_cupy.view(cp.float16)
mask = af_float16 > cp.float16(0.5)
start_run = time.perf_counter()
df22 = df2[mask.get()]
df22 = df22[df22["AC"] > 8]
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & AF > 0.5 & GT = 1|1 & AD > 3: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

#Filter AC > 8 & AF > 0.5 & GT = 1|1 & AD[0] > 3 & AD[1] < 10
af_cupy = cp.asarray(df2["AF"].values)
af_float16 = af_cupy.view(cp.float16)
mask = af_float16 > cp.float16(0.5)
start_run = time.perf_counter()
df22 = df2[mask.get()]
df22 = df22[df22["AC"] > 8]
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3) & (df3["AD1"] < 10)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Tempo  AC > 8 & AF > 0.5 & GT = 1|1 & AD[0] > 3 & AD[1] < 10: {time_elapsed:.4f} secondi")
tot_time = tot_time + time_elapsed

del df33
del df22
gc.collect()

print(f"Tempo totale di query: {tot_time:.4f} secondi")