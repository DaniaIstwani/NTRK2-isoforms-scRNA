

# 📄 GTF Transcript Extractor

This script extracts specific **transcript entries** from a GTF file using transcript IDs. It's useful for isolating annotations related to specific genes or transcripts — for example, extracting isoforms of a gene like `Rbms1`.

---

## 🧬 Purpose

Given a list of transcript IDs, this script:
- Extracts all lines in the GTF file that mention those transcript IDs.
- Optionally includes surrounding lines for context (all parent gene entries).
- Removes duplicate lines to produce a clean, minimal file.

---

# 🔍 Extract Rbms1 Transcript Entries from GTF

This script extracts specific **Rbms1 transcript entries** (e.g., `Rbms1-001`, `Rbms1-209`) from a GTF annotation file. It includes both the **transcript lines** and their **associated parent gene lines**, if present.

---

## 🧾 Script: `extract_rbms1_transcripts.sh`

```bash
#!/bin/bash

# Input and output files
input_gtf="original.gtf"  # Your original GTF file
output_txt="rbms1_transcripts.txt"  # Output file with extracted lines

# Transcript IDs (adjust if needed)
rbms1_001="ENSMUST00000028347.12"
rbms1_209="ENSMUST00000164147.7"

# Extract all lines containing these transcript IDs
grep -E "transcript_id \"($rbms1_001|$rbms1_209)\"" "$input_gtf" > "$output_txt"

# Also include their parent 'gene' lines (optional)
grep -A 1 -B 1 -E "transcript_id \"($rbms1_001|$rbms1_209)\"" "$input_gtf" >> "$output_txt"

# Remove duplicate lines (if any)
sort -u "$output_txt" > temp && mv temp "$output_txt"

echo "Extracted all lines for Rbms1-001 and Rbms1-209 to: $output_txt"
```



# 🛠️ Transcript Attribute Cleaner and Modifier

This script processes a subset of GTF transcript entries — for example, those extracted for `Rbms1` — by simplifying and modifying the attributes column. It rewrites transcript entries to represent minimal custom gene models such as `short-Rbms1` and `long-Rbms1`.

---

## 🎯 Purpose

Given a GTF file containing selected transcript entries (e.g., from `rbms1_transcripts.txt`), this script:

1. Replaces the `gene_id` field with the **transcript ID**.
2. Assigns a custom **gene name** (`short-Rbms1` or `long-Rbms1`) based on the transcript ID.
3. Removes the `transcript_id` and **all other attributes**, retaining only:
   - `gene_id`
   - `gene_name`
4. Converts feature type `"transcript"` to `"gene"` (for simplification).

---

## 🧾 Script: `process_rbms1_transcripts.sh`

```bash
#!/bin/bash

input_file="rbms1_transcripts.txt"
output_file="processed_rbms1.gtf"

# Process the file
awk -F'\t' '
BEGIN {OFS = FS}
{
    # Extract transcript_id
    tid = ""
    split($9, attrs, ";")
    for (i in attrs) {
        if (attrs[i] ~ /transcript_id/) {
            split(attrs[i], t, "\"")
            tid = t[2]
            break
        }
    }

    # Skip if no transcript_id found
    if (tid == "") {
        print
        next
    }

    # Set gene name based on transcript_id
    gname = (tid == "ENSMUST00000164147.7") ? "short-Rbms1" : "long-Rbms1"

    # Rebuild attributes
    new_attrs = "gene_id \"" tid "\"; gene_name \"" gname "\";"

    # Update fields
    $9 = new_attrs
    if ($3 == "transcript") {
        $3 = "gene"
    }

    print
}' "$input_file" > "$output_file"

echo "Processing complete. Output saved to $output_file"
```

# 🧬 Final GTF Cleaner & Custom Transcript Merger

This script removes all traces of a specific gene (e.g., `Rbms1`) from a GTF annotation file and inserts custom-curated transcript entries in its place. The result is a clean, final GTF annotation ready for visualization or downstream analysis.

---

## 🎯 Purpose

Given:
- An **original GTF** file (e.g., from Ensembl or GENCODE)
- A **custom GTF** containing processed transcripts (e.g., `short-Rbms1`, `long-Rbms1`)

This script:
1. Removes all entries associated with the original gene and its transcripts.
2. Appends curated transcript entries from a custom GTF.
3. Sorts the result to maintain GTF compatibility.
4. Validates the replacement.
5. Outputs a clean `final_annotation.gtf`.

---

## 🧾 Script: `finalize_gtf.sh`

```bash
#!/bin/bash

# File paths
original_gtf="original.gtf"
processed_transcripts="processed_rbms1.gtf"
final_gtf="final_annotation.gtf"

# === 1. Remove all traces of Rbms1 ===
grep -v -E \
    -e 'gene_id "ENSMUSG00000026970"' \
    -e 'gene_name "Rbms1"' \
    -e 'transcript_name "Rbms1' \
    -e 'transcript_id "ENSMUST00000164147\.7"' \
    -e 'transcript_id "ENSMUST00000028347\.12"' \
    "$original_gtf" > temp_stage1.gtf

# === 2. Secondary sweep to remove case-insensitive "rbms1" remnants ===
grep -v -i "rbms1" temp_stage1.gtf > temp_stage2.gtf

# === 3. Add custom-processed Rbms1 transcripts ===
cat "$processed_transcripts" >> temp_stage2.gtf

# === 4. Sort entries, preserving any GTF header lines ===
{
    grep "^#" temp_stage2.gtf 2>/dev/null || true
    grep -v "^#" temp_stage2.gtf | sort -k1,1 -k4,4n
} > "$final_gtf"

# === 5. Validation ===
echo "=== VALIDATION ==="
echo "Original Rbms1 count: $(grep -ic "rbms1" "$original_gtf")"
echo "Final Rbms1 count (excluding custom names): $(grep -i "rbms1" "$final_gtf" | grep -viE 'short|long' | wc -l)"
echo "Custom entries added:"
echo "  short-Rbms1: $(grep -c 'gene_name \"short-Rbms1\"' "$final_gtf")"
echo "  long-Rbms1:  $(grep -c 'gene_name \"long-Rbms1\"' "$final_gtf")"

# === 6. Cleanup ===
rm temp_stage*.gtf
echo "✅ Process completed. Final file: $final_gtf"
```
