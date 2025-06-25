#!/usr/bin/env python3

#source cyvcf2_env/bin/activate


from cyvcf2 import VCF, Writer
import time

# Definisci funzioni filtro per ogni condizione
def filter_EVA_4(variant):
    # Restituisce True se il campo INFO "EVA_4" esiste (non None)
    return variant.INFO.get("EVA_4") is not None

def filter_E_Multiple_observations(variant):
    return variant.INFO.get("E_Multiple_observations") is not None

def filter_TSA_SNV(variant):
    # Supponiamo che il campo "TSA" contenga la stringa "SNV"
    tsa = variant.INFO.get("TSA")
    return tsa is not None and tsa.strip() == "SNV"

def filter_pos_gt(variant, pos_min):
    return variant.POS > pos_min

def filter_pos_range(variant, pos_min, pos_max):
    return pos_min < variant.POS < pos_max

def filter_not_EMultiple_and_TSA_SNV(variant):
    """
    Passa la variante se:
      - Il campo INFO "E_Multiple_observations" NON è presente (None),
      - E il campo "TSA" esiste e, pulito dagli spazi, è uguale a "SNV".
    """
    has_EMultiple = variant.INFO.get("E_Multiple_observations")
    tsa = variant.INFO.get("TSA")
    return (has_EMultiple is None) and (tsa is not None and tsa.strip() == "SNV")

def filter_pos_gt_and_EMultiple(variant):
    """
    Passa la variante se (su file annotato):
      - Il campo INFO "mypos" è presente ed è maggiore di 1780000,
      - E il campo INFO "E_Multiple_observations" è presente.
    """
    mypos = variant.INFO.get("mypos")
    has_EMultiple = variant.INFO.get("E_Multiple_observations")
    try:
        pos_val = int(mypos) if mypos is not None else None
    except ValueError:
        pos_val = None
    return (pos_val is not None and pos_val > 1780000) and (has_EMultiple is not None)

def filter_pos_range_and_EMultiple(variant):
    """
    Passa la variante se (su file annotato):
      - Il campo INFO "mypos" è presente ed è compreso tra 1780000 e 1800000,
      - E il campo INFO "E_Multiple_observations" è presente.
    """
    mypos = variant.INFO.get("mypos")
    has_EMultiple = variant.INFO.get("E_Multiple_observations")
    try:
        pos_val = int(mypos) if mypos is not None else None
    except ValueError:
        pos_val = None
    return (pos_val is not None and 1780000 < pos_val < 1800000) and (has_EMultiple is not None)

def filter_range_and_TSA_SNV(variant):
    """
    Filtro: POS > 1780000 & POS < 1800000 & TSA=SNV - danio_rerio
    Si applica su un file annotato contenente il campo INFO "mypos".
    Restituisce True se:
      - "mypos" esiste ed il suo valore (convertito a int) è compreso tra 1780000 e 1800000,
      - E il campo INFO "TSA" esiste e, rimosso eventuale spazi, è uguale a "SNV".
    """
    mypos = variant.INFO.get("mypos")
    tsa = variant.INFO.get("TSA")
    try:
        pos_val = int(mypos) if mypos is not None else None
    except ValueError:
        pos_val = None
    return (pos_val is not None and 1780000 < pos_val < 1800000) and (tsa is not None and tsa.strip() == "SNV")

def filter_EVA4_range_and_TSA_SNV(variant):
    """
    Filtro: EVA_4 & POS > 1780000 & POS < 1800000 & TSA=SNV - danio_rerio
    Restituisce True se:
      - Il campo INFO "EVA_4" esiste,
      - "mypos" esiste ed il suo valore (convertito a int) è compreso tra 1780000 e 1800000,
      - E il campo INFO "TSA" esiste e, rimosso eventuale spazi, è uguale a "SNV".
    """
    eva4 = variant.INFO.get("EVA_4")
    mypos = variant.INFO.get("mypos")
    tsa = variant.INFO.get("TSA")
    try:
        pos_val = int(mypos) if mypos is not None else None
    except ValueError:
        pos_val = None
    return (eva4 is not None) and (pos_val is not None and 1780000 < pos_val < 1800000) and (tsa is not None and tsa.strip() == "SNV")


# Funzione generale per applicare un filtro a un VCF
def filter_vcf(input_vcf, output_vcf, filter_func):
    vcf_reader = VCF(input_vcf)
    vcf_writer = Writer(output_vcf, vcf_reader)
    for variant in vcf_reader:
        if filter_func(variant):
            vcf_writer.write_record(variant)
    vcf_writer.close()

def main():
    # Esempio filtro: EVA_4
    print("Filtro: EVA_4")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_EVA_4.vcf", filter_EVA_4)
    end = time.perf_counter()
    print("Tempo EVA_4: {:.4f} s".format(end - start))
    
    # Filtro: E_Multiple_observations
    print("Filtro: E_Multiple_observations")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_E_Multiple_observations.vcf", filter_E_Multiple_observations)
    end = time.perf_counter()
    print("Tempo E_Multiple_observations: {:.4f} s".format(end - start))
    
    # Filtro: TSA = SNV
    print("Filtro: TSA=SNV")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_TSA_SNV.vcf", filter_TSA_SNV)
    end = time.perf_counter()
    print("Tempo TSA=SNV: {:.4f} s".format(end - start))
    
    # Filtro: POS > 1780000
    print("Filtro: POS > 1780000")
    start = time.perf_counter()
    # Qui usiamo il campo POS direttamente (non è nel INFO, perciò non serve annotazione)
    filter_vcf("data/danio_rerio.vcf", "output_POS_gt_1780000.vcf",
               lambda var: filter_pos_gt(var, 1780000))
    end = time.perf_counter()
    print("Tempo POS > 1780000: {:.4f} s".format(end - start))
    
    # Filtro: POS tra 1780000 e 1800000
    print("Filtro: 1780000 < POS < 1800000")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_POS_1780000_1800000.vcf",
               lambda var: filter_pos_range(var, 1780000, 1800000))
    end = time.perf_counter()
    print("Tempo POS range: {:.4f} s".format(end - start))

    # Filtro: ! E_Multiple_observations & TSA = SNV (da VCF originale)
    print("Filtro: ! E_Multiple_observations & TSA = SNV")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_not_EMultiple_TSA_SNV.vcf", filter_not_EMultiple_and_TSA_SNV)
    end = time.perf_counter()
    print("Tempo ! E_Multiple_observations & TSA = SNV: {:.4f} s".format(end - start))
    
    # Filtro: POS > 1780000 & E_Multiple_observations (su file annotato)
    print("Filtro: POS > 1780000 & E_Multiple_observations")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio_annotated_new_fixed.vcf", "output_POS_gt_and_EMultiple.vcf",
               filter_pos_gt_and_EMultiple)
    end = time.perf_counter()
    print("Tempo POS > 1780000 & E_Multiple_observations: {:.4f} s".format(end - start))

    # Filtro: POS > 1780000 & POS < 1800000 & E_Multiple_observations (su file annotato)
    print("Filtro: 1780000 < POS < 1800000 & E_Multiple_observations")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio_annotated_new_fixed.vcf", "output_POS_range_and_EMultiple.vcf",
               filter_pos_range_and_EMultiple)
    end = time.perf_counter()
    print("Tempo 1780000 < POS < 1800000 & E_Multiple_observations: {:.4f} s".format(end - start))

    # Filtro combinato: POS > 1780000 & TSA=SNV
    print("Filtro: POS > 1780000 & TSA=SNV")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_combined.vcf",
               lambda var: var.POS > 1780000 and filter_TSA_SNV(var))
    end = time.perf_counter()
    print("Tempo  POS > 1780000 & TSA=SNV: {:.4f} s".format(end - start))

    # Filtro: POS > 1780000 & POS < 1800000 & TSA=SNV - danio_rerio
    print("Filtro: 1780000 < POS < 1800000 & TSA = SNV (danio_rerio)")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio_annotated_new_fixed.vcf", "output_range_TSA_SNV_bos.vcf",
               filter_range_and_TSA_SNV)
    end = time.perf_counter()
    print("Tempo 1780000 < POS < 1800000 & TSA = SNV: {:.4f} s".format(end - start))
        
    # Filtro combinato: EVA_4 & POS > 1780000 & TSA=SNV
    print("Filtro: EVA_4 & POS > 1780000 & TSA=SNV")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio.vcf", "output_combined.vcf",
               lambda var: filter_EVA_4(var) and var.POS > 1780000 and filter_TSA_SNV(var))
    end = time.perf_counter()
    print("Tempo filtro combinato: {:.4f} s".format(end - start))

    # Filtro: EVA_4 & POS > 1780000 & POS < 1800000 & TSA=SNV - danio_rerio
    print("Filtro: EVA_4 & 1780000 < POS < 1800000 & TSA = SNV (danio_rerio)")
    start = time.perf_counter()
    filter_vcf("data/danio_rerio_annotated_new_fixed.vcf", "output_EVA4_range_TSA_SNV_bos.vcf",
               filter_EVA4_range_and_TSA_SNV)
    end = time.perf_counter()
    print("Tempo EVA_4 & 1780000 < POS < 1800000 & TSA = SNV: {:.4f} s".format(end - start))
    
    
if __name__ == "__main__":
    main()
