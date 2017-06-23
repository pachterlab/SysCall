


## learns coefficients from coverage file
## runs on list
## outputs to one file sys errors
## outputs to second file hetero sites. 

## to work with NaN values in PT. Doesn't train on that. 

rm(list = ls())
library(MASS);

if (!(7==length(commandArgs()))){
        cat("ERROR: not right number of variables \n");
        stop("USAGE: variables train_file out_sys out_hetero \n");
}

instances <-commandArgs()[4];
train_file <-commandArgs()[5];
outf_sys <-commandArgs()[6];
outf_hetero <-commandArgs()[7];



train.df = read.table(train_file, header=T)

train.df$IS_SYS_ERROR = factor(train.df$IS_SYS_ERROR, labels=c(1,0))

attach(train.df)

logit.1 = glm(IS_SYS_ERROR  ~ base.2 +base.1 + base + diff_err_diff_dir + diff_error_dir , family=binomial(link=logit))


test.df = read.table(instances, header=T)

pred = predict.glm(logit.1, newdata=test.df, type="response");

pred<-1-pred;
bind_pred<-cbind(test.df, pred);

sys_errs<-na.omit(bind_pred[bind_pred$pred<=0.5, c(1,length(bind_pred[1,]))]);

hetero_locations<-na.omit(bind_pred[bind_pred$pred>0.5, c(1,length(bind_pred[1,]))]);

if ((length(hetero_locations[,1])+length(sys_errs[,1]))!= length(bind_pred[,1])){
	stop("ERROR in lengths \n");	
}

write.table(sys_errs,file=outf_sys,quote=F,row.names=F,col.names=F);
write.table(hetero_locations,file= outf_hetero,quote=F,row.names=F,col.names=F);

detach(train.df)



