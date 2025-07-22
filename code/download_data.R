library(misha.ext)
gdb.create_genome("mm10", here("data"))
gdb.create_genome("hg38", here("data"))

download.file("https://iceqream.s3.eu-west-1.amazonaws.com/gastrulation_energies.tar.gz", "data/gastrulation_energies.tar.gz")
untar("data/gastrulation_energies.tar.gz")