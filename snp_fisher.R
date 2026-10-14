d <- read.delim("D:/claude-study/Spinal Ependymoma/R_pkg_logs/_snp_overlap.tsv", header=FALSE,
                col.names=c("name","listn","inmat","ndmp","ndmp_total","n1000","a","b","c","d"))
for (i in 1:nrow(d)) {
  t2 <- matrix(c(d$a[i], d$b[i], d$c[i], d$d[i]), 2)
  ft <- fisher.test(t2)
  cat(d$name[i], sprintf("OR %.3f p %.4g\n", ft$estimate, ft$p.value))
}
