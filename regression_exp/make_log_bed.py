import math

in_bed="/home/asm242/EntexExpress/entex_data/ENCFF646XTO.bed"
out_bed="peaks.bed"

max_seq_len = -1

with open(in_bed) as infile, open(out_bed, "w") as outfile:
    for line in infile:
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = fields[1]
        end = fields[2]

        if abs(int(end)-int(start)) > max_seq_len or max_seq_len == 0:
            max_seq_len = abs(int(end)-int(start))

        signal = float(fields[6])
        log_signal = math.log(signal + 1e-6)  # small epsilon to avoid log(0)

        name = f"|chr={chrom}|signal={signal:.3f}|log(signal)={log_signal:.3f}|"
        out_fields = [chrom, start, end, name, f"{signal:.5f}", f"{log_signal:.5f}"]
        outfile.write("\t".join(out_fields) + "\n")

print(max_seq_len)

