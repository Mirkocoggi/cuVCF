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

# Parsing e misurazione dei tempi
res = vcf.vcf_parsed()
print("Start parsing")
res.run("data/IRBT.vcf", 16)

data1 = vcf.get_var_columns_data(res.var_columns) 
data2 = vcf.get_alt_columns_data(res.alt_columns)
data3 = vcf.get_sample_columns_data(res.samp_columns)
data4 = vcf.get_alt_format_data(res.alt_sample)

# Converti "qual" da uint16 a float32, se necessario
#if "qual" in data1:
#    arr = data1["qual"].view(np.float16)
#    data1["qual"] = arr.astype(np.float32)

# parch perchè non va il binding su var_id
n = len(data3["var_id"])
group_size = 8
data3["var_id"] = np.repeat(np.arange((n + group_size - 1) // group_size), group_size)[:n]

df1 = cudf.DataFrame(data1)
df2 = cudf.DataFrame(data2)
df3 = cudf.DataFrame(data3)
df4 = cudf.DataFrame(data4)

print("Dataframe creati")

filter_map = vcf.GTMapGlobal
ref_value = filter_map["1|1"]

#Filter AC>8
start_run = time.perf_counter()
df22 = df2[df2["AC"] > 8]
end_run = time.perf_counter()
print(f"Filter AC>8: {end_run - start_run:.4f} secondi")
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
print(f"Filter AC > 8 & AF > 0.5: {end_run - start_run:.4f} secondi")
print(len(df2))
print(len(df22))
del df22
gc.collect()

#Filter GT = 1|1
start_run = time.perf_counter()
df33 = df3[df3["GT0"] == ref_value]
end_run = time.perf_counter()
print(f"Filter GT = 1|1: {end_run - start_run:.4f} secondi")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AD[0] > 3
start_run = time.perf_counter()
df33 = df3[df3["AD0"] > 3]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3: {end_run - start_run:.4f} secondi")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AD[0] > 3 & AD[1] < 10
start_run = time.perf_counter()
df33 = df3[(df3["AD0"] > 3) & (df3["AD1"] < 10)]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & AD[1] < 10: {end_run - start_run:.4f} secondi")
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
print(f"Filter GT = 1|1 & AC>8: {end_run - start_run:.4f} secondi")
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
print(f"Filter GT = 1|1 & AC>8 & AF > 0.5: {end_run - start_run:.4f} secondi")
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
print(f"Filter AD[0] > 3 & AC>8: {end_run - start_run:.4f} secondi")
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
print(f"Filter AC > 8 & AD[0] > 3 & AD[1] < 10: {end_run - start_run:.4f} secondi")
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
print(f"Filter AC > 8 & AF > 0.5 & AD[0] > 3: {end_run - start_run:.4f} secondi")
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
print(f"Filter AC > 8 & AF > 0.5 & AD[0] > 3 & AD[1] < 10: {end_run - start_run:.4f} secondi")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()

#Filter AD[0] > 3 & GT = 1|1
start_run = time.perf_counter()
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3)]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & GT = 1|1: {end_run - start_run:.4f} secondi")
print(len(df3))
print(len(df33))
del df33
gc.collect()

#Filter AD[0] > 3 & AD[1] < 10 & GT = 1|1
start_run = time.perf_counter()
df33 = df3[(df3["GT0"] == ref_value) & (df3["AD0"] > 3) & (df3["AD1"] < 10)]
end_run = time.perf_counter()
print(f"Filter AD[0] > 3 & AD[1] < 10 & GT = 1|1: {end_run - start_run:.4f} secondi")
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
print(f"Filter AC > 8 & GT = 1|1 & AD[0] > 3: {end_run - start_run:.4f} secondi")
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
print(f"Filter AC > 8 & AD[0] > 3 & AD[1] < 10 & GT = 1|1: {end_run - start_run:.4f} secondi")
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
print(f"Filter AAC > 8 & AF > 0.5 & GT = 1|1 & AD > 3: {end_run - start_run:.4f} secondi")
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
print(f"Filter AAC > 8 & AF > 0.5 & GT = 1|1 & AD > 3: {end_run - start_run:.4f} secondi")
print(len(df3))
print(len(df33))
del df33
del df22
gc.collect()


#def filter_complete_genotypes(df3, threshold=0.75, total_samples=8):
#    """
#    Filtra df3 eliminando le varianti in cui la percentuale di genotipi validi,
#    ossia quelli il cui valore in "GT0" non è il carattere corrispondente a missing (chr(255)),
#    è inferiore al 'threshold'. Con total_samples=8, viene mantenuta la variante solo se almeno
#    ceil(8 * threshold) chiamate sono valide.
#    
#    Parameters:
#        df3 (cudf.DataFrame): DataFrame con almeno le colonne "var_id" e "GT0".
#        threshold (float): Percentuale minima di completezza, default 0.75.
#        total_samples (int): Numero totale di campioni per variante, default 8.
#    
#    Returns:
#        cudf.DataFrame: DataFrame filtrato con solo le varianti per cui il numero di genotipi validi
#                        è >= ceil(total_samples * threshold).
#    """
#    # Calcola il numero minimo di genotipi validi richiesto (es. ceil(8*0.75)=6)
#    min_valid = int(np.ceil(total_samples * threshold))
#    
#    # Il codice missing è rappresentato da chr(255)
#    missing_code = chr(255)
#    
#    # Crea una colonna booleana "valid": True se "GT0" è diverso dal codice missing
#    df_temp = df3.assign(valid=(df3["GT0"] != missing_code))
#    
#    # Raggruppa per "var_id" e conta quanti "valid" (True vale 1) sono presenti per ciascuna variante
#    valid_counts = df_temp.groupby("var_id").agg({"valid": "sum"}).reset_index()
#    # Seleziona i var_id che hanno almeno min_valid genotipi validi
#    valid_ids = valid_counts[valid_counts["valid"] >= min_valid]["var_id"]
#    
#    # Filtra df3 mantenendo solo le varianti con var_id in valid_ids
#    return df3[df3["var_id"].isin(valid_ids)]
#
## Filtra df3 escludendo le varianti con percentuale di genotipi validi < 0.75
#start_run = time.perf_counter()
#df3_filtered = filter_complete_genotypes(df3, 0.75, 8)
#end_run = time.perf_counter()
#print(f"Filter genotipi validi: {end_run - start_run:.4f} secondi")
#print(len(df3))
#print(len(df3_filtered))
#del df3_filtered
#gc.collect()

#def filter_alt_depth(df3, depth_threshold=8, min_ratio=0.5):
#    """
#    Filtra il DataFrame df3 (che contiene le colonne "var_id", "AD1", ecc.)
#    escludendo le varianti in cui la percentuale di campioni con profondità
#    dell’allele alternativo (AD1) maggiore o uguale a depth_threshold risulti inferiore a min_ratio.
#    
#    Parameters:
#        df3 (cudf.DataFrame): DataFrame dei genotipi, contenente almeno le colonne "var_id" e "AD1".
#        depth_threshold (int): Soglia minima per AD1 (default 8).
#        min_ratio (float): Percentuale minima di campioni (default 0.5) che devono superare la soglia.
#        
#    Returns:
#        cudf.DataFrame: Nuovo DataFrame filtrato contenente solo le varianti per cui
#                        la condizione è verificata.
#                        
#    Nota: df3 non viene modificato.
#    """
#    # Crea un DataFrame temporaneo con una colonna booleana "valid"
#    # che indica se AD1 del campione è >= depth_threshold.
#    df_temp = df3.assign(valid=(df3["AD1"] >= depth_threshold))
#    
#    # Raggruppa per "var_id" per ottenere, per ogni variante, il numero di
#    # chiamate valide e il numero totale di campioni.
#    grouped = df_temp.groupby("var_id").agg({"valid": "sum", "AD1": "count"}).reset_index()
#    grouped = grouped.rename(columns={"valid": "valid_count", "AD1": "total_count"})
#    
#    # Filtra le varianti in cui la proporzione di chiamate valide è >= min_ratio
#    valid_ids = grouped[grouped["valid_count"] / grouped["total_count"] >= min_ratio]["var_id"]
#    
#    # Restituisci un nuovo DataFrame filtrato, senza modificare df3
#    return df3[df3["var_id"].isin(valid_ids)]

# Filtra df3 escludendo le varianti per con percentuale di campioni AD1 >= 8 inferiore a 0.5
#start_run = time.perf_counter()
#df3_filtered = filter_alt_depth(df3, 8, 0.5)
#end_run = time.perf_counter()
#print(f"Filter AD1 >= 8 nel 50% samp: {end_run - start_run:.4f} secondi")
#print(len(df3))
#print(len(df3_filtered))
#del df3_filtered
#gc.collect()
#
#def filter_window_position(df1, chrom_target=0, pos_min=150, pos_max=3000):
#    """
#    Filtra il DataFrame df1 per mantenere solo le varianti localizzate sul cromosoma "1"
#    (codificato come 0) e con la posizione compresa nell'intervallo [pos_min, pos_max].
#
#    Parameters:
#        df1 (cudf.DataFrame): DataFrame contenente i dati varianti con almeno le colonne "chrom" e "pos".
#        chrom_target (int): Valore target per identificare il cromosoma "1". Default è 0.
#        pos_min (int): Posizione minima (default 150).
#        pos_max (int): Posizione massima (default 3000).
#
#    Returns:
#        cudf.DataFrame: Un nuovo DataFrame contenente solo le varianti che soddisfano i criteri di filtro.
#    """
#    return df1[(df1["chrom"] == chrom_target) & 
#               (df1["pos"] >= pos_min) & 
#               (df1["pos"] <= pos_max)]

# Filtro per l'intervallo genomico (cromosoma "1" e posizione [150, 3000])
#start_run = time.perf_counter()
#df1_filtered = filter_window_position(df1)
#end_run = time.perf_counter()
#print(f"Filter cromosoma 1 e posizione [150, 3000]: {end_run - start_run:.4f} secondi")
#print(len(df1))
#print(len(df1_filtered))
#del df1_filtered
#gc.collect()


#filter 

#print("Prime 3 righe di df1:")
#print(df1.head(3))
#
#print("Prime 3 righe di df2:")
#print(df2.head(3))
#
#print("Prime 3 righe di df3:")
#print(df3.head(3))
#
#print("Prime 3 righe di df4:")
#print(df4.head(3))

#npartitions = 10  # puoi regolare questo parametro in base alle tue necessità
#save_cudf_to_csv_in_chunks(df1, "df1.csv", npartitions, index=False)
#save_cudf_to_csv_in_chunks(df2, "df2.csv", npartitions, index=False)
#save_cudf_to_csv_in_chunks(df3, "df3.csv", npartitions, index=False)
#save_cudf_to_csv_in_chunks(df4, "df4.csv", npartitions, index=False)

