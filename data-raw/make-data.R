count <- read.table("data-raw/combo_count_tab.txt")
depth <- sapply(strsplit(colnames(count), "\\."), length)
x <- count[, depth == 6 & colSums(count != 0) >= 1]

demo <- read.delim("data-raw/demographic.txt")
y <- demo$bmi[match(rownames(count), demo$pid)]
bmi_gut = data.frame(bmi = y, x)

save(bmi_gut, file = "data/bmi_gut.rda", compress = "bzip2")
