## Two proportion test of conservation percentage between overlapping groups and non-overlapping groups.

## Parallel overlap
prop.test(x=c(12,436), n=c(22,996), 
          alternative = c("two.sided"),
          conf.level = 0.95, correct = TRUE)

## Convergentl overlap
prop.test(x=c(274,436), n=c(371,996), 
          alternative = c("two.sided"),
          conf.level = 0.95, correct = TRUE)

## Divergent overlap -> sample size less than 5 need to use fisher's exact test
testor=rbind(c(2,436),c(13,996))
fisher.test(testor,alternative=c("two.sided"))


## Parent-embed anti-parallel overlap
prop.test(x=c(55,436), n=c(94,996), 
          alternative = c("two.sided"),
          conf.level = 0.95, correct = TRUE)

## Parent-embed parallel overlap
prop.test(x=c(14,436), n=c(69,996), 
          alternative = c("two.sided"),
          conf.level = 0.95, correct = TRUE)
