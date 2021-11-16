file1 <- "/home/jfsanchez/DATA/XICRA/test_read_len/seqs_test"
file2 <- "/home/jfsanchez/DATA/XICRA/test_read_len/all_seqs.freqs.subset_isomiRs.fasta"

?readDNAStringSet

## dataframe
# index, seq, len, name, biotype, group, count

## migth be able to plot according to colour extracted from header

############################
## readDNAStringSet
############################
## fasta file or fastq format file
library(Biostrings)
file2_seqs <- readDNAStringSet(filepath = file2, format = "fasta", use.names = TRUE)

seqs_df <- as.data.frame(x = file2_seqs, row.names = NULL)
seqs_df$name <- row.names(seqs_df)
rownames(seqs_df) <- NULL
colnames(seqs_df)[1] <- "seqs"
seqs_df$len <- nchar(seqs_df$seqs)

library(dplyr)
seqs_df_filtered <- seqs_df %>% 
  group_by(seqs) %>% 
  mutate(IDs = paste0(unique(name), collapse = ",")) %>% count(IDs, name = "count") %>% as.data.frame()

seqs_df_filtered$len <- nchar(seqs_df_filtered$seqs)
seqs_df_filtered$biotype <- "miRNA"
seqs_df_filtered$group <- "sample1"

head(seqs_df_filtered)
dim(seqs_df_filtered)
############################

############################
## just sequences, no fasta header
############################
file1_seqs <- readLines(file1)

name_list <- as.character()
for(i in 1:length(file1_seqs)){
  name_list <- append(name_list, paste0('seq_',i))
}
name_list

seqs_df2 <- as.data.frame(x = file1_seqs, row.names = NULL)
colnames(seqs_df2)[1] <- "seqs"
seqs_df2$len <- nchar(seqs_df2$seqs)
seqs_df2$name <- name_list

seqs_df_filtered2 <- seqs_df2 %>% 
  group_by(seqs) %>% 
  mutate(IDs = paste0(unique(name), collapse = ",")) %>% count(IDs, name = "count") %>% as.data.frame()

seqs_df_filtered2$len <- nchar(seqs_df_filtered2$seqs)
seqs_df_filtered2$biotype <- "miRNA"
seqs_df_filtered2$group <- "sample2"

seqs_df_filtered2 <- seqs_df_filtered2[seqs_df_filtered2$len<=25,]

dim(seqs_df_filtered2)
dim(seqs_df_filtered)

############################

############################
## plot
############################
# base R
hist(seqs_df$len)

# ggplot
library(ggplot2)
ggplot(seqs_df, mapping = aes(len)) + geom_bar()
ggplot(seqs_df_filtered, mapping = aes(x=len)) + geom_bar() ## unique

ggplot(seqs_df_filtered, mapping = aes(x=len, y=count)) + geom_bar(stat="identity") + ## all
  labs(title = "Frequency by isomiR length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic()


## example biotypes
seqs_df3_filtered <- seqs_df %>% 
  group_by(seqs) %>% 
  mutate(IDs = paste0(unique(name), collapse = ",")) %>% count(IDs, name = "count") %>% as.data.frame()

head(seqs_df3_filtered)
seqs_df3_filtered$len <- nchar(seqs_df3_filtered$seqs)

seqs_df3_filtered$biotype <- "miRNA"
seqs_df3_filtered$group <- ""

seqs_df4_filtered <- seqs_df3_filtered[sample(nrow(seqs_df3_filtered), 400),]
seqs_df4_filtered$biotype <- "tRNA"
seqs_df4_filtered$len <- seqs_df4_filtered$len*4
seqs_df4_filtered$count <- seqs_df4_filtered$count*2

dim(seqs_df3_filtered)
dim(seqs_df4_filtered)

seqs_df5_filtered <- rbind(seqs_df3_filtered, seqs_df4_filtered)
ggplot(seqs_df5_filtered, mapping = aes(x=len, y=count, fill=biotype)) + geom_bar(stat="identity") + ## all
  labs(title = "Frequency by sequence length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() #+ facet_grid(biotype ~ .)


# example different samples
seqs_df6_filtered <- rbind(seqs_df_filtered, seqs_df_filtered2)
ggplot(seqs_df6_filtered, mapping = aes(x=len, y=count, fill=group)) + geom_bar(stat="identity") + ## all
  labs(title = "Frequency by sequence length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() + facet_grid(group ~ .)

ggplot(seqs_df6_filtered, mapping = aes(x=len, y=count, fill=group)) + 
  geom_bar(stat="identity") + ## all
  labs(title = "Frequency by sequence length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() #+ facet_grid(group ~ .)


ggplot(seqs_df6_filtered, mapping = aes(x=len, y=count, fill=group)) + geom_bar(stat="identity") + ## all
  labs(title = "Frequency by sequence length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() #+ facet_grid(group ~ .)


ggplot(seqs_df6_filtered, mapping = aes(x=len, y=count, fill=group, colour=group)) + 
  stat_ecdf() +
  geom_bar(stat="identity") + ## all
  labs(title = "Frequency by sequence length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() #+ facet_grid(group ~ .)


  


p1 <- seqs_df6_filtered %>>% ggplot() +  
  geom_bar(mapping = aes(x=len, y=count, fill=group), stat="identity") + 
  stat_ecdf(mapping = aes(x=len, y=count, colour=group), 
               sec.axis = sec_axis(~ . * 1/4000, name = "Relative humidity [%]"))
  
p1

labs(title = "Frequency by sequence length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() #+ facet_grid(group ~ .)

