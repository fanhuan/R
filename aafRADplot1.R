## Part of package Huan
## aafRADplot1.R

# Read in SBA data
#RAD<-read.delim("~/Dropbox/Research/RAD/Condor_results/result_160613_RAD_d01.txt", col.names=c('ID','L','r','ks','mis_b','mis_sba','opt_k','mis_a','kmers_b','kmers_sba','kmers_a'))
# Read in collapsed pairwise data
#RAD_d01_pairwise<- read.delim("~/Dropbox/Research/RAD/Condor_results/160613_RAD_d01_pairwise_ks.txt", header=FALSE, col.names = c('ID', 'L', 'r', 'ks', 'k', 'mis_p'))
#Collapse and merge data

aafRADplot1 <- function (df_SBA,df_pairwise,title) {
    total_col<-ddply(df_SBA,~ID+L+r,summarise,mis_b=min(mis_b),mis_a=min(mis_a),mis_sba=min(mis_sba))
    total_col_1 <- merge(total_col,df_pairwise, by = c('ID','L','r'))
    total_col_mean<-total_col_1 %>% 
        group_by(L,r) %>%
        summarize(se_b=sd(mis_b)/sqrt(length(mis_b)),se_a=sd(mis_a)/sqrt(length(mis_a)),se_sba=sd(mis_sba)/sqrt(length(mis_sba)),se_p=sd(mis_p)/sqrt(length(mis_p)),mis_b=mean(mis_b),mis_a=mean(mis_a),mis_sba=mean(mis_sba),mis_p=mean(mis_p)) 
    limit_a <- aes(ymax = mis_a + se_a, ymin=mis_a - se_a)
    limit_sba <- aes(ymax = mis_sba + se_sba, ymin=mis_sba - se_sba)
    limit_b <- aes(ymax = mis_b + se_b, ymin=mis_b - se_b)
    limit_p <- aes(ymax = mis_p + se_p, ymin=mis_p - se_p)
    #Plot
    ggplot(total_col_mean,aes(x=L,y=mis_b)) +
        ylab('topological mistakes') +
        xlab('number of loci') +
        ggtitle(title) +
        facet_wrap(~r) +
        geom_point(pch='b',size=4) +
        geom_errorbar(limit_b, width=1) +
        geom_line() +
        geom_point(aes(y=mis_a),pch='a',size=4,color='blue') +
        geom_errorbar(limit_a, width=0.2,color='blue') +
        geom_line(aes(y=mis_a),color='blue',lty=2) +
        geom_point(aes(y=mis_sba),pch='s',size=4,color='red') +
        geom_errorbar(limit_sba, width=0.2,color='red') +
        geom_line(aes(y=mis_sba),color='red',lty=3) +
        geom_point(aes(y=mis_p),pch='p',size=4,color='green') +
        geom_errorbar(limit_p, width=0.2,color='green') +
        geom_line(aes(y=mis_p),color='green',lty=5) 
}
