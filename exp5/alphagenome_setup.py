from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import pandas as pd

# Load gtf from hg38 gencode
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

# Set up transcript extractors
gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)

# Access API with personal API key
dna_model = dna_client.create("AIzaSyDaeU_HLHvVAN6tpNK3D82dTVx1olTzxuI")

# Pick 2048 bp around ENSEMBL gene ID ENSG00000188976
interval = gene_annotation.get_gene_interval(gtf, gene_id='ENSG00000188976')
interval = interval.resize(dna_client.SEQUENCE_LENGTH_2KB)

# Predict RNA_SEQ in Brain Substantia Nigra in this interval
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0002038'], # Brain_Substantia_Nigra
)

# Extract longest transcripts & plot them
longest_transcripts = transcript_extractor.extract(interval)
print(f'Extracted {len(longest_transcripts)} transcripts in this interval.')

plot_components.plot(
    components=[
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.Tracks(output.rna_seq),
    ],
    interval=output.rna_seq.interval,
)

plt.savefig("brain_sn_rna_seq.png")
plt.close()



