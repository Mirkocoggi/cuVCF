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

res = vcf.vcf_parsed()
print("Start parsing")
res.run("data/IRBT.vcf", 16)

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

print("Dataframes created")

filter_map = vcf.GTMapGlobal
ref_value = filter_map["1|1"]

#Filter AC>8
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
end_run = time.perf_counter()
print(f"Filter AC>8: {end_run - start_run:.4f} seconds")
print(len(df2))
print(len(df22))
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
print(f"Filter AC > 8 & AF > 0.5: {end_run - start_run:.4f} seconds")
print(len(df2))
print(len(df22))
del df22
gc.collect()

#Filter GT = 1|1
start_run = time.perf_counter()
df33 = df3[df3["GT0"] == ref_value]
end_run = time.perf_counter()
print(f"Filter GT = 1|1: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AD[0] > 3
start_run = time.perf_counter()
df33 = df3[df3["AD0"] > 3]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AD[0] > 3 & AD[1] < 10
start_run = time.perf_counter()
df33 = df3[(df3["AD0"] > 3) & (df3["AD1"] < 10)]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & AD[1] < 10: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter GT = 1|1 & AC>8
start_run = time.perf_counter()
df33 = df3[df3["GT0"] == ref_value]
df22 = df2[df2["AC"] > 8]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
print(f"Filter GT = 1|1 & AC>8: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
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
print(f"Filter GT = 1|1 & AC>8 & AF > 0.5: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()

#Filter AD[0] > 3 & AC>8
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[df3["AD0"] > 3]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & AC>8: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()

#Filter AC > 8 & AD[0] > 3 & AD[1] < 10
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[(df3["AD0"] > 3) & (df3["AD1"] < 10)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
print(f"Filter AC > 8 & AD[0] > 3 & AD[1] < 10: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
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
print(f"Filter AC > 8 & AF > 0.5 & AD[0] > 3: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
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
print(f"Filter AC > 8 & AF > 0.5 & AD[0] > 3 & AD[1] < 10: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()

#Filter AD[0] > 3 & GT = 1|1
start_run = time.perf_counter()
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3)]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & GT = 1|1: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AD[0] > 3 & AD[1] < 10 & GT = 1|1
start_run = time.perf_counter()
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3) & (df3["AD1"] < 10)]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & AD[1] < 10 & GT = 1|1: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AC > 8 & GT = 1|1 & AD[0] > 3
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
print(f"Filter AC > 8 & GT = 1|1 & AD[0] > 3: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()

#Filter AC > 8 & AD[0] > 3 & AD[1] < 10 & GT = 1|1
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3) & (df3["AD1"] < 10)]
df33 = df33[df33["var_id"].isin(df22["var_id"])]
end_run = time.perf_counter()
print(f"Filter AC > 8 & AD[0] > 3 & AD[1] < 10 & GT = 1|1: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
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
print(f"Filter AAC > 8 & AF > 0.5 & GT = 1|1 & AD > 3: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
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
print(f"Filter AAC > 8 & AF > 0.5 & GT = 1|1 & AD > 3: {end_run - start_run:.4f} seconds")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()
