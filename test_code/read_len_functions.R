sessionInfo()


file1 <- "./seqs_test"
file2 <- "./seqs_test.fasta"

?readDNAStringSet

## dataframe
# index, len, biotype, group, count, [name, seq]

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
## just data in tab format
file3 <- "/imppc/labs/lslab/jsanchez/DATA/XICRA/test_RNAbiotype/test_sample.length.results.txt"
data3 <- read.table(file3, sep="\t")
data3

############################
## plot
############################
# base R
hist(seqs_df$len)

# ggplot

install.packages("ggplot2")

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
############################


############################
## just data in tab format
file3 <- "/imppc/labs/lslab/jsanchez/DATA/XICRA/test_RNAbiotype/test_sample.length.results.txt"
data3 <- read.table(file3, sep="\t", header = TRUE)
dim(data3)

table(data3$biotype)
hist(data3$length)
data4 <- head(data3)

# ggplot
library(ggplot2, lib.loc = "/soft/general/R-3.6.3-bioc-3.10/lib/R/library/")

colnames(data3)

ggplot(data3, mapping = aes(x=length)) + geom_bar() ## unique

ggplot(data3, mapping = aes(x=length, y=count)) + geom_bar(stat="identity") + ## all
  labs(title = "Frequency by isomiR length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic()

ggplot(data3, mapping = aes(x=length, y=count, fill=biotype)) + geom_bar(stat="identity") + ## all
  labs(title = "Frequency by isomiR length\n", x = "\nLength (bp)", y = "Frequency\n") +
  theme_classic() + facet_grid(biotype ~ .)

# ---------------------------------------------------------------------------------
## summarize
# ---------------------------------------------------------------------------------

biotypes <- data.frame(type=data3$biotype,count=data3$count)
sum.biotypes<- aggregate(biotypes[,2], by=list(biotypes$type), "sum")
sum.biotypes$Group.1 <- factor(sum.biotypes$Group.1, levels=unique(biotypes$type))

sum.biotypes

## generate percentage
perc.biotypes<- cbind.data.frame(sum.biotypes$Group.1, 100 * prop.table( as.matrix(sum.biotypes[, 2:ncol(sum.biotypes)]), 2 ) )
names(perc.biotypes)[1]<- "Biotypes"
names(perc.biotypes)[2]<- "percentage"

library(dplyr)
major_categories <- perc.biotypes %>%
  filter(percentage >= 5) %>%
  pull("Biotypes")

data3 %>%
  filter(biotype %in% major_categories) %>%
  ggplot(mapping = aes(x=length, y=count, fill=biotype)) + geom_bar(stat="identity") + ## all
    labs(title = "Frequency by Biotype class\n", x = "\nLength (bp)", y = "Frequency\n") +
    theme_classic() + facet_grid(biotype ~ .)



