print("starting script")
library(arrow)

args <- commandArgs(trailingOnly = TRUE)


input_files <- head(args, -1)
print("Input files: ")
print(input_files)
output_file = tail(args, 1)
print(paste0("Output file: ", output_file))

print("Reading chunked rds files and combining them")
for (file in input_files) {

    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")) {
        dataset <- readRDS(file)
    }

    # if the merged dataset does exist, append to it
    if (exists("dataset")) {
        temp_dataset <- readRDS(file)
        dataset <- rbind(dataset, temp_dataset)
        rm(temp_dataset)
    }
}


dataset["gene"] = sapply(unlist(dataset["gene"]), as.numeric)

print(paste0("Writing compbined ouptut to ", output_file))
write_parquet(dataset, sink = paste0(output_file))

