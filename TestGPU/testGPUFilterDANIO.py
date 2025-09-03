# source rapids_env/bin/activate
import cupy as cp
import GPUParser as vcf
import cudf
import time
import numpy as np
import math
import gc

def save_cudf_to_csv_in_chunks(df, filename, npartitions=10, index=False):

    n_rows = len(df)
    chunk_size = math.ceil(n_rows / npartitions)
    
    header = ",".join(df.columns.astype(str)) + "\n"
    
    with open(filename, "w") as f:
        f.write(header)
        for i in range(npartitions):
            start = i * chunk_size
            end = min(start + chunk_size, n_rows)
            chunk = df.iloc[start:end]
            csv_str = chunk.to_csv(index=index, header=False)
            f.write(csv_str)

time_elapsed = 0

res = vcf.vcf_parsed()
print("Start parsing")
start_run = time.perf_counter()
res.run("data/danio_rerio.vcf", 16)
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Parsing time: {time_elapsed:.4f} seconds")

start_run = time.perf_counter()

data1 = vcf.get_var_columns_data(res.var_columns) 
data2 = vcf.get_alt_columns_data(res.alt_columns)
data3 = vcf.get_sample_columns_data(res.samp_columns)
data4 = vcf.get_alt_format_data(res.alt_sample)

n = len(data3["var_id"])
group_size = 8
data3["var_id"] = np.repeat(np.arange((n + group_size - 1) // group_size), group_size)[:n]
df1 = cudf.DataFrame(data1)
df2 = cudf.DataFrame(data2)
df3 = cudf.DataFrame(data3)
df4 = cudf.DataFrame(data4)
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"Dataframes created in: {time_elapsed:.4f} seconds")

start_run = time.perf_counter()
npartitions = 10  
save_cudf_to_csv_in_chunks(df1, "df1.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df2, "df2.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df3, "df3.csv", npartitions, index=False)
save_cudf_to_csv_in_chunks(df4, "df4.csv", npartitions, index=False)
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"CSV created in: {time_elapsed:.4f} seconds")

del df1
del df2
del df3
del df4

start_run = time.perf_counter()
df1 = cudf.read_csv("df1.csv", delimiter=",")
df2 = cudf.read_csv("df2.csv", delimiter=",")
df3 = cudf.read_csv("df3.csv", delimiter=",")
df4 = cudf.read_csv("df4.csv", delimiter=",")
end_run = time.perf_counter()
time_elapsed = time_elapsed+(end_run - start_run)
print(f"CSV loaded in: {time_elapsed:.4f} seconds")

#Filter EVA_4
start_run = time.perf_counter()
df11 = df1[df1["EVA_4"] == 1]
end_run = time.perf_counter()
print(f"Filter EVA_4: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter E_Multiple_observations
start_run = time.perf_counter()
df11 = df1[df1["E_Multiple_observations"] == 1]
end_run = time.perf_counter()
print(f"Filter E_Multiple_observations: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter TSA = SNV
start_run = time.perf_counter()
df11 = df1[df1["TSA"] == 0]
end_run = time.perf_counter()
print(f"Filter TSA = SNV: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter POS > 1780000
start_run = time.perf_counter()
df11 = df1[df1["pos"] > 1780000]
end_run = time.perf_counter()
print(f"Filter POS > 1780000: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter POS > 1780000 & POS < 1800000
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["pos"] < 1800000)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & POS < 1800000: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter ! E_Multiple_observations and TSA=SNV 
start_run = time.perf_counter()
df11 = df1[(df1["E_Multiple_observations"] == 0) & (df1["TSA"] == 0)]
end_run = time.perf_counter()
print(f"Filter ! E_Multiple_observations and TSA = SNV: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter POS > 1780000 & E_Multiple_observations 
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["E_Multiple_observations"] == 1)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & E_Multiple_observations: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter POS > 1780000 & POS < 1800000 & E_Multiple_observations
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["pos"] < 1800000) & (df1["E_Multiple_observations"] == 1)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & POS < 1800000 & E_Multiple_observations: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter POS > 1780000 & TSA=SNV 
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["TSA"] == 0)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & TSA=SNV : {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter POS > 1780000 & POS < 1800000 & TSA=SNV
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["pos"] < 1800000) & (df1["TSA"] == 0)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & POS < 1800000 & TSA=SNV: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter EVA_4 & POS > 1780000 & TSA=SNV
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["TSA"] == 0) & (df1["EVA_4"] == 1)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & TSA=SNV : {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()

#Filter EVA_4 & POS > 1780000 & POS < 1800000 & TSA=SNV
start_run = time.perf_counter()
df11 = df1[(df1["pos"] > 1780000) & (df1["pos"] < 1800000) & (df1["TSA"] == 0) & (df1["EVA_4"] == 1)]
end_run = time.perf_counter()
print(f"Filter POS > 1780000 & POS < 1800000 & TSA=SNV: {end_run - start_run:.4f} seconds")
print(len(df1))
print(len(df11))
del df11
gc.collect()