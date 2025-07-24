library(misha.ext)

options(timeout = 900)

# Download and extract the core IQ paper data first
download.file("https://iq-paper.s3.ap-south-1.amazonaws.com/IQ-paper-data.tar.gz", "IQ-paper-data.tar.gz")
untar("IQ-paper-data.tar.gz", exdir = ".")

# Create genome databases
gdb.create_genome("mm10", here("data"))
gdb.create_genome("hg38", here("data"))

# Download and extract tracks for hg38
download.file("https://iq-paper.s3.ap-south-1.amazonaws.com/hg38-tracks.tar.gz", "data/hg38-tracks.tar.gz")
untar("data/hg38-tracks.tar.gz", exdir = "data/hg38/")

# Download and extract tracks for mm10
download.file("https://iq-paper.s3.ap-south-1.amazonaws.com/mm10-tracks.tar.gz", "data/mm10-tracks.tar.gz")
untar("data/mm10-tracks.tar.gz", exdir = "data/mm10/")

# Add the gastrulation energies data
download.file("https://iceqream.s3.eu-west-1.amazonaws.com/gastrulation_energies.tar.gz", "data/gastrulation_energies.tar.gz")
untar("data/gastrulation_energies.tar.gz")
