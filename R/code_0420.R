##check if MB better than optimal Bon
##take alpha as input, compute sample size (should be at least smaller than opt bon)
##binary search for n still start at simple bon

#result: opt bon vs. MB vs. MB using input alpha
#when p=0.4, 232 vs. 250 vs.231
#when p=0.7, 298 vs. 310 vs. 299

IfPowerSatisfy.fix.alpha <- function(pi.1, delta1, delta2, sigma, N, beta, alpha.input){
        ###functions
        #mean and covariate matrix under alternative
        ConvertToMuSigma <- function(pi.1,delta1,delta2,sigma,N){
                pi.1 <- pi.1
                delta.1 <- delta1
                delta.2 <- delta2
                sigma <- sigma
                #compute Sigma
                pi.2 <- 1-pi.1
                Sigma <- matrix(c(1,sqrt(pi.1),sqrt(pi.2),sqrt(pi.1),1,0,sqrt(pi.2),0,1),nrow = 3,ncol = 3)
                row.names(Sigma) <- c(0,1,2)
                colnames(Sigma) <- c(0,1,2)
                #compute mu
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                N.0 <- N #total sample size
                N.1 <- N.0*pi.1
                N.2 <- N.0*pi.2
                #alternative under 3 different power condition
                #(,1,1)
                mu.0 <- c(delta.0/(2*sigma/sqrt(N.0)),delta.1/(2*sigma/sqrt(N.1)),delta.2/(2*sigma/sqrt(N.2)))
                #(,1,0)
                delta.1 <- delta1
                delta.2 <- 0
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                mu.1 <- c(delta.0/(2*sigma/sqrt(N.0)),delta.1/(2*sigma/sqrt(N.1)),delta.2/(2*sigma/sqrt(N.2)))
                #(,0,1)
                delta.1 <- 0
                delta.2 <- delta2
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                mu.2 <- c(delta.0/(2*sigma/sqrt(N.0)),delta.1/(2*sigma/sqrt(N.1)),delta.2/(2*sigma/sqrt(N.2)))
                return(list(mu = cbind(mu.0,mu.1,mu.2), Sigma = Sigma))
        }
        #convert to C
        alpha.matrix <- function(alpha.star){
                C <- matrix(0,nrow = 7,ncol = 3)
                type1.error <- sum(alpha.star)
                C[1,1] = C[2,2] = C[3,3] = type1.error
                #test for (c,1)
                C[4,1] = alpha.star[1]
                C[4,2] = alpha.star[2]+alpha.star[3]
                #test for (c,2)
                C[5,1] = alpha.star[1]
                C[5,3] = alpha.star[2]+alpha.star[3]
                #test for (1,2)
                C[6,2] = alpha.star[2]+1/2*alpha.star[1]
                C[6,3] = alpha.star[3]+1/2*alpha.star[1]
                #test for (c,1,2)
                C[7,1] = alpha.star[1]
                C[7,2] = alpha.star[2]
                C[7,3] = alpha.star[3]
                rownames(C) <- c("0","1","2","0,1","0,2","1,2","0,1,2")
                colnames(C) <- c("0","1","2")
                return(C)
        }
        #calculate power
        PowerMB <- function(C,mu,Sigma,j){
                C <- C
                mu <- mu #a matrix under different alternative
                Sigma <- Sigma
                ####Function: local rejection region for H_F
                RejRegionHFBon <- function(F){
                        for (i in 1:nrow(C)){
                                if(row.names(C)[i] == F){
                                        c_j <- C[i,]
                                        break;
                                }
                        }
                        alpha.current <- c_j
                        u <- numeric(3)
                        index <- as.numeric(unlist(strsplit(F,",")))
                        index <- index+1
                        for (i in 1:length(index)){
                                u[index[i]] <- qnorm(alpha.current[index[i]],mean = 0, sd = 1,lower.tail = FALSE)
                        }
                        return(u)
                }

                ####Function: rejection region for rejecting a single test by closed testing procedure, i.e reject all local tests          ####containing this single test
                RejRegionSingle <- function(j){
                        cut.point <- NULL #store local rejection region, 4 rows
                        for (i in 1:nrow(C)){#local rejection regions
                                if(j %in% unlist(strsplit(rownames(C)[i],","))){
                                        F <- rownames(C)[i]
                                        cut.point[[i]] <- RejRegionHFBon(F)
                                }
                        }
                        #combine all rejection region together
                        cut.point <- data.frame(t(matrix(unlist(cut.point),nrow = 3))) #equivalent with A vector!
                        colnames(cut.point) <- c("0","1","2")

                        ####Function:Calculate all rejection cubes of intersection of unions
                        RejCube <- function(cut.point){
                                index <- as.numeric(j)+1
                                index.others <- setdiff(c(1,2,3),index)
                                represent.point <- NULL
                                #construct for test 1
                                a <- cut.point[,index]
                                a <- a[order(a)]
                                temp <- NULL
                                temp[1] <- a[1]-1
                                temp[2] <- (a[1]+a[2])/2
                                temp[3] <- (a[2]+a[3])/2
                                temp[4] <- (a[3]+a[4])/2
                                temp[5] <- a[4]+1
                                represent.point[[index]] <- temp
                                #construct for test 0
                                index <- index.others[1]
                                a <- cut.point[,index][cut.point[,index]!=0]
                                a <- a[order(a)]
                                temp <- NULL
                                temp[1] <- a[1]-1
                                temp[2] <- (a[1]+a[2])/2
                                temp[3] <- a[2]+1
                                represent.point[[index]] <- temp
                                #constrct for test 2
                                index <- index.others[2]
                                a <- cut.point[,index][cut.point[,index]!=0]
                                a <- a[order(a)]
                                temp <- NULL
                                temp[1] <- a[1]-1
                                temp[2] <- (a[1]+a[2])/2
                                temp[3] <- a[2]+1
                                represent.point[[index]] <- temp

                                IfPointInLocalRegion <- function(p,i){ #i is the row index in cut.point
                                        index <- which(cut.point[i,]!=0)
                                        if(any(p[index]>cut.point[i,][index]))#OR
                                                return(TRUE)
                                        else return(FALSE)
                                }

                                IfPointInRegion <- function(p){
                                        temp <- 1
                                        for (i in 1:nrow(cut.point)){
                                                temp <- temp * IfPointInLocalRegion(p,i)
                                        }
                                        if (temp==1) return(TRUE)
                                        else return(FALSE)
                                }

                                cube <- NULL
                                x <- cut.point[,1][cut.point[,1]!=0]
                                x <- x[order(x)]
                                x <- append(append(x,Inf),-Inf,after=0)
                                y <- cut.point[,2][cut.point[,2]!=0]
                                y <- append(append(y,Inf),-Inf,after=0)
                                y <- y[order(y)]
                                z <- cut.point[,3][cut.point[,3]!=0]
                                z <- z[order(z)]
                                z <- append(append(z,Inf),-Inf,after=0)

                                id <- 1
                                for(i in 1:length(represent.point[[1]])){#for loop
                                        for(j in 1:length(represent.point[[2]])){
                                                for(k in 1:length(represent.point[[3]])){
                                                        p <- c(represent.point[[1]][i],represent.point[[2]][j],represent.point[[3]][k])
                                                        if(IfPointInRegion(p)){
                                                                cube[[id]] <- rbind(c(x[i],x[i+1]),c(y[j],y[j+1]),c(z[k],z[k+1]))
                                                                id <- id+1
                                                        }
                                                }
                                        }
                                }#end of for
                                return(cube)
                        }
                        return(RejCube(cut.point))
                }

                ####Function: compute the exact power for rejecting test j under alternative mu_j
                PowerRejSingle <- function(j,mu,Sigma){
                        shape <- RejRegionSingle(j) #rejection region
                        number.of.cube <- length(shape)
                        area <- 0
                        for (i in 1:number.of.cube){
                                lower <- shape[[i]][,1]
                                upper <- shape[[i]][,2]
                                area <- area + pmvnorm(lower = lower, upper = upper, mean = mu[,(as.numeric(j)+1)], sigma = Sigma,algorithm=GenzBretz(abseps = 10^-10 ,maxpts=10^5))
                        }
                        return(as.numeric(area))
                }

                return(PowerRejSingle(j,mu,Sigma))
        }

        mu <- ConvertToMuSigma(pi.1,delta1,delta2,sigma,N)$mu
        Sigma <- ConvertToMuSigma(pi.1,delta1,delta2,sigma,N)$Sigma

        C <- alpha.matrix(alpha.input)
        if((PowerMB(C,mu,Sigma,"0") > beta) & (PowerMB(C,mu,Sigma,"1") > beta) & (PowerMB(C,mu,Sigma,"2") > beta)) return(list(satisfy=TRUE))
        else return(list(satisfy=FALSE))
}


SampleSizeMB.fix.alpha <- function(pi.1,delta.1,delta.2,sigma,alpha,beta,alpha.input){
        SampleSizeBon <- function(pi.1,delta.1,delta.2,sigma,alpha,beta){
                alpha <- alpha/2
                alpha <- alpha/3 #bon adjustment, one-sided test
                pi.1 <- pi.1
                pi.2 <- 1-pi.1
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                delta.1 <- delta.1
                delta.2 <- delta.2
                sigma <- sigma
                n0.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.0^2)
                n1.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.1^2)
                n2.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.2^2)
                #satisfy all 3 power conditions
                return(ceiling(max(n0.req,n1.req/pi.1,n2.req/pi.2)))
        }
        #search from size.bon
        n.upper <- floor(SampleSizeBon(pi.1,delta.1,delta.2,sigma,alpha,beta))
        n.lower <- 0
        n.current <- ceiling((n.upper +n.lower)/2)
        #current mu and Sigma
        repeat{
                temp <- IfPowerSatisfy.fix.alpha(pi.1,delta.1,delta.2,sigma,n.current,beta,alpha.input)
                if(temp$satisfy){
                        n.upper <- n.current
                        n.current <- ceiling((n.upper +n.lower)/2)
                        alpha.current <- temp$alpha
                }
                if(!temp$satisfy){
                        n.lower <- n.current
                        n.current <- ceiling((n.upper +n.lower)/2)
                }
                if(n.upper-n.lower <= 1) break
        }

        return(list(sample.size=n.upper))
}

###################simulation

###FUNCTION: calculate rejection region
RejRegion <- function(alpha.star,j){
        ###functions
        #convert to C
        alpha.matrix <- function(alpha.star){
                C <- matrix(0,nrow = 7,ncol = 3)
                type1.error <- sum(alpha.star)
                C[1,1] = C[2,2] = C[3,3] = type1.error
                #test for (c,1)
                C[4,1] = alpha.star[1]
                C[4,2] = alpha.star[2]+alpha.star[3]
                #test for (c,2)
                C[5,1] = alpha.star[1]
                C[5,3] = alpha.star[2]+alpha.star[3]
                #test for (1,2)
                C[6,2] = alpha.star[2]+1/2*alpha.star[1]
                C[6,3] = alpha.star[3]+1/2*alpha.star[1]
                #test for (c,1,2)
                C[7,1] = alpha.star[1]
                C[7,2] = alpha.star[2]
                C[7,3] = alpha.star[3]
                rownames(C) <- c("0","1","2","0,1","0,2","1,2","0,1,2")
                colnames(C) <- c("0","1","2")
                return(C)
        }
        #calculate rejection region
        RejRegionSingle <- function(C,j){
                C <- C
                RejRegionHFBon <- function(F){
                        for (i in 1:nrow(C)){
                                if(row.names(C)[i] == F){
                                        c_j <- C[i,]
                                        break;
                                }
                        }
                        alpha.current <- c_j
                        u <- numeric(3)
                        index <- as.numeric(unlist(strsplit(F,",")))
                        index <- index+1
                        for (i in 1:length(index)){
                                u[index[i]] <- qnorm(alpha.current[index[i]],mean = 0, sd = 1,lower.tail = FALSE)
                        }
                        return(u)
                }
                cut.point <- NULL #store local rejection region, 4 rows
                for (i in 1:nrow(C)){#local rejection regions
                        if(j %in% unlist(strsplit(rownames(C)[i],","))){
                                F <- rownames(C)[i]
                                cut.point[[i]] <- RejRegionHFBon(F)
                        }
                }
                #combine all rejection region together
                cut.point <- data.frame(t(matrix(unlist(cut.point),nrow = 3))) #equivalent with A vector!
                colnames(cut.point) <- c("0","1","2")

                ####Function:Calculate all rejection cubes of intersection of unions
                RejCube <- function(cut.point){
                        index <- as.numeric(j)+1
                        index.others <- setdiff(c(1,2,3),index)
                        represent.point <- NULL
                        #construct for test 1
                        a <- cut.point[,index]
                        a <- a[order(a)]
                        temp <- NULL
                        temp[1] <- a[1]-1
                        temp[2] <- (a[1]+a[2])/2
                        temp[3] <- (a[2]+a[3])/2
                        temp[4] <- (a[3]+a[4])/2
                        temp[5] <- a[4]+1
                        represent.point[[index]] <- temp
                        #construct for test 0
                        index <- index.others[1]
                        a <- cut.point[,index][cut.point[,index]!=0]
                        a <- a[order(a)]
                        temp <- NULL
                        temp[1] <- a[1]-1
                        temp[2] <- (a[1]+a[2])/2
                        temp[3] <- a[2]+1
                        represent.point[[index]] <- temp
                        #constrct for test 2
                        index <- index.others[2]
                        a <- cut.point[,index][cut.point[,index]!=0]
                        a <- a[order(a)]
                        temp <- NULL
                        temp[1] <- a[1]-1
                        temp[2] <- (a[1]+a[2])/2
                        temp[3] <- a[2]+1
                        represent.point[[index]] <- temp

                        IfPointInLocalRegion <- function(p,i){ #i is the row index in cut.point
                                index <- which(cut.point[i,]!=0)
                                if(any(p[index]>cut.point[i,][index]))#OR
                                        return(TRUE)
                                else return(FALSE)
                        }

                        IfPointInRegion <- function(p){
                                temp <- 1
                                for (i in 1:nrow(cut.point)){
                                        temp <- temp * IfPointInLocalRegion(p,i)
                                }
                                if (temp==1) return(TRUE)
                                else return(FALSE)
                        }

                        cube <- NULL
                        x <- cut.point[,1][cut.point[,1]!=0]
                        x <- x[order(x)]
                        x <- append(append(x,Inf),-Inf,after=0)
                        y <- cut.point[,2][cut.point[,2]!=0]
                        y <- append(append(y,Inf),-Inf,after=0)
                        y <- y[order(y)]
                        z <- cut.point[,3][cut.point[,3]!=0]
                        z <- z[order(z)]
                        z <- append(append(z,Inf),-Inf,after=0)

                        id <- 1
                        for(i in 1:length(represent.point[[1]])){#for loop
                                for(j in 1:length(represent.point[[2]])){
                                        for(k in 1:length(represent.point[[3]])){
                                                p <- c(represent.point[[1]][i],represent.point[[2]][j],represent.point[[3]][k])
                                                if(IfPointInRegion(p)){
                                                        cube[[id]] <- rbind(c(x[i],x[i+1]),c(y[j],y[j+1]),c(z[k],z[k+1]))
                                                        id <- id+1
                                                }
                                        }
                                }
                        }#end of for
                        return(cube)
                }
                return(RejCube(cut.point))
        }

        C <- alpha.matrix(alpha.star)
        rej.region <- RejRegionSingle(C,j)
        rej.region
}

###FUNCTION: if z statistics fall in the rejection region
IfInRejRegion <- function(rej.region,z){
        IfinCube <- function(i){
                cube <- rej.region[[i]]
                cube1 <-cube[1,]
                cube2 <-cube[2,]
                cube3 <-cube[3,]
                if((z[1]>=cube1[1])&(z[1]<=cube1[2])&(z[2]>=cube2[1])&(z[2]<=cube2[2])&(z[3]>=cube3[1])&(z[3]<=cube3[2])) return(TRUE)
                else return(FALSE)
        } #judge if in the single cube
        for (i in 1:length(rej.region)){
                if (IfinCube(i)) return(TRUE) #if in one of the single cube, return
                else next() #if not, judge for next cube
        }
        #if not in any of the cubes, return false
        return(FALSE)
}

