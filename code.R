library(abind)
library(tidyverse)
library(doParallel)

### Function ####

deom_ediff = function(yeard,l0_1,l0_2,P_1,P_2){
  
  D = dim(l0_1)[1]
  
  lx <- function(trans, init){
    l = list(init)
    for (a in 1:(length(trans))) {
      lx = l[[a]] %*% trans[[a]]
      l[[length(l)+1]]=lx
    }
    return(l)
  }
  
  l_1 = lx(P_1,l0_1)
  l_2 = lx(P_2,l0_2)
  
  HLE <- function(trans, init){
    l = list(init)
    e = init*0
    for (a in 1:(length(trans))) {
      lx = l[[a]] %*% trans[[a]]
      l[[length(l)+1]]=lx
      Lx = (l[[a]]+l[[a+1]])/2
      e = e+Lx
    }
    return(e)
  }
  
  e_1 = HLE(P_1,l0_1)
  e_2 = HLE(P_2,l0_2)
  
  
  m_l = list()
  d_l = list()
  for (i in 1:length(l_1)) {
    m_l[[i]] = (l_1[[i]]+l_2[[i]])/2
    d_l[[i]] = (l_2[[i]]-l_1[[i]])/yeard
  }
  
  m_P = list()
  d_P = list()
  for (i in 1:length(P_1)) {
    m_P[[i]] = (P_1[[i]]+P_2[[i]])/2
    d_P[[i]] = (P_2[[i]]-P_1[[i]])/yeard
  }
  
  
  m_e = (e_1+e_2)/2
  d_e = (e_2-e_1)/yeard
  
  #
  All_MoP = vector(mode = "list",length = length(m_P))
  for (i in 1:(length(m_P)-1)) {
    All_MoP[[i]] = list(m_P[[i]])
    for (j in i:(length(m_P)-1)) {
      tem = All_MoP[[i]][[length(All_MoP[[i]])]] %*% m_P[[j+1]]
      All_MoP[[i]][[length(All_MoP[[i]])+1]] = tem
    }
    All_MoP[[i]][[length(All_MoP[[i]])]] = All_MoP[[i]][[length(All_MoP[[i]])]]/2
  }
  All_MoP[[(length(m_P))]] = list(m_P[[(length(m_P))]]/2)
  
  I= diag(D)
  
  # all Ps after the p-dot
  SUM = list()
  for(i in 1:length(All_MoP)){
    SUM[[i]]= Reduce(`+`, All_MoP[[i]]) + I/2
  }
  for(i in 2:length(SUM)){
    SUM[[i]]= SUM[[i]] + I/2
  }
  SUM[[length(SUM)+1]]=I/2
  
  # lx*px_dot
  Internal = vector(mode = "list",length = length(SUM))
  DOT = list(d_l[[1]])
  for (i in 2:length(SUM)) {
    DOT[[i]] = m_l[[i-1]] %*% d_P[[i-1]]
    # for further dientangle in the Px_dot
    internal = c()
    for(j in c(1:D)){
      for (k in c(1:D)) {
        for (g in c(1:D)) {
          tem = m_l[[i-1]][j,]*d_P[[i-1]][,g]*SUM[[i]][g,k]
          internal = append(internal,tem)
        }
      }
    }
    
    # reshape the data and give col & row names
    internal = matrix(internal,D,D^3,byrow = T)
    cn = c()
    for (j in c(1:D)) {
        for (g in c(1:D)) {
          cn = append(cn, paste0(g,j))
        }
    }
    A = array(0,c(D,D,D^2),dimnames = list(c(1:D),c(1:D),cn))
    for (j in c(1:D^2)) {
      A[,,j] = internal[,c(seq(j,D^3,by=D^2))]
    }
    
    Internal[[i]] = A
    
  }
  Internal[[1]] = array(0,c(D,D,D^2),dimnames = list(c(1:D),c(1:D),cn))
  
  res = Map('%*%',DOT,SUM)
  
  RES = list(Input=list(P_1=P_1,P_2=P_2,l0_1=l0_1,l0_2=l0_2,yeard=yeard),e_1=e_1,e_2=e_2,d_e=d_e,res=res,d_p = Internal)
  RES
}



### Data ####
INI = read_csv(paste0("RMLE/BASELINE.csv"))
TRANS = read_csv(paste0("RMLE/PROB.csv"))

# Input ####
yeard = 1 # t2-t1 (time difference)

registerDoParallel(min(detectCores(),12))
ALL_RES = foreach(i=1:500,.packages = c("tidyverse")) %dopar% {
  
  ## l0 in t1 and t2
  ini = INI %>% filter(ragender == 1,iter==i)%>% ungroup()%>% 
    arrange(state) %>% pull(pro)
  
  t= (ini[1]+ini[2])
  ini[1] = ini[1]/t
  ini[2] = ini[2]/t
  l0_1 = matrix(c(ini[1],0,0,
                  0,ini[2],0,
                  0,0,0),nrow=3,byrow = T)
  
  ini = INI %>% filter(ragender == "2",iter==i)%>% ungroup()%>% 
    arrange(state) %>% pull(pro)
  
  t= (ini[1]+ini[2])
  ini[1] = ini[1]/t
  ini[2] = ini[2]/t
  l0_2 = matrix(c(ini[1],0,0,
                  0,ini[2],0,
                  0,0,0),nrow=3,byrow = T)
  
  ## transition probabilities of men
  trans = TRANS %>% filter(ragender == 1,iter==i) %>% 
    pivot_longer(c(5:7),names_to = "state",values_to = "prob")
  H = crossing(pre_state= "H",state = c("A","L","H"),ragender= unique(trans$ragender),
               age = unique(trans$age),iter = i,prob = 0)
  H = H %>% mutate(prob = ifelse(state == "H"& pre_state=="H", 1, prob))
  trans = bind_rows(trans,H)
  trans$pre_state = factor(trans$pre_state,levels = c("A","L","H"))
  trans$state = factor(trans$state,levels = c("A","L","H"))
  trans = xtabs(prob ~ pre_state+state+age,data = trans)

  P_1 = lapply(seq(dim(trans)[3]), function(x) trans[ , , x])
  
  ## transition probabilities of women
  trans = TRANS %>% filter(ragender == "2",iter==i) %>% 
    pivot_longer(c(5:7),names_to = "state",values_to = "prob")
  H = crossing(pre_state= "H",state = c("A","L","H"),ragender= unique(trans$ragender),
               age = unique(trans$age),iter = i,prob = 0)
  H = H %>% mutate(prob = ifelse(state == "H"& pre_state=="H", 1, prob))
  trans = bind_rows(trans,H)
  trans$pre_state = factor(trans$pre_state,levels = c("A","L","H"))
  trans$state = factor(trans$state,levels = c("A","L","H"))
  trans = xtabs(prob ~ pre_state+state+age,data = trans)
  
  P_2 = lapply(seq(dim(trans)[3]), function(x) trans[ , , x])
  
  P_1 = lapply(P_1, function(X){
    X<- X/rowSums(X)
    X
  })
  
  P_2 = lapply(P_2, function(X){
    X<- X/rowSums(X)
    X
  })
  
  ##make it 2x2
  # cutdeath = function(x){
  #   x = x[1:2,1:2]
  # }
  # 
  # l0_1 = cutdeath(l0_1)
  # l0_2 = cutdeath(l0_2)
  # 
  # P_1 = lapply(P_1, cutdeath)
  # P_2 = lapply(P_2, cutdeath)
  
  res = deom_ediff(yeard,l0_1,l0_2,P_1,P_2)
  res
}
stopImplicitCluster()

### Results ####
e_1 = ALL_RES[[1]][["e_1"]]
for (a in 2:500) {
  e_1 = abind(e_1,ALL_RES[[a]][["e_1"]],along = 3)
}
round(ALL_RES[[1]][["e_1"]],2)
round(colSums(ALL_RES[[1]][["e_1"]]),2)
round(apply(e_1,c(1,2),quantile,0.025),2)
round(apply(e_1,c(1,2),quantile,0.975),2)
round(apply(apply(e_1,c(2,3),sum),c(1),quantile,0.025),2)
round(apply(apply(e_1,c(2,3),sum),c(1),quantile,0.975),2)

e_2 = ALL_RES[[1]][["e_2"]]
for (a in 2:500) {
  e_2 = abind(e_2,ALL_RES[[a]][["e_2"]],along = 3)
}
round(ALL_RES[[1]][["e_2"]],2)
round(colSums(ALL_RES[[1]][["e_2"]]),2)
round(apply(e_2,c(1,2),quantile,0.025),2)
round(apply(e_2,c(1,2),quantile,0.975),2)
round(apply(apply(e_2,c(2,3),sum),c(1),quantile,0.025),2)
round(apply(apply(e_2,c(2,3),sum),c(1),quantile,0.975),2)


e_dot = ALL_RES[[1]][["d_e"]]
for (a in 2:500) {
  e_dot = abind(e_dot,ALL_RES[[a]][["d_e"]],along = 3)
}
round(ALL_RES[[1]][["d_e"]],2)
round(colSums(ALL_RES[[1]][["d_e"]]),3)
round(apply(e_dot,c(1,2),quantile,0.025),3)
round(apply(e_dot,c(1,2),quantile,0.975),3)
round(apply(apply(e_dot,c(2,3),sum),c(1),quantile,0.025),3)
round(apply(apply(e_dot,c(2,3),sum),c(1),quantile,0.975),3)


l_dot = ALL_RES[[1]][["res"]][[1]]
for (a in 2:500) {
  l_dot = abind(l_dot,ALL_RES[[a]][["res"]][[1]],along = 3)
}
round(ALL_RES[[1]][["res"]][[1]],2)
round(colSums(ALL_RES[[1]][["res"]][[1]]),2)
round(apply(apply(l_dot,c(2,3),sum),c(1),quantile,0.025),3)
round(apply(l_dot,c(1,2),quantile,0.025),3)
# round(apply(l_dot,c(1,2),quantile,0.5),3)
round(apply(l_dot,c(1,2),quantile,0.975),3)
round(apply(apply(l_dot,c(2,3),sum),c(1),quantile,0.975),3)

P_dot = Reduce(`+`, ALL_RES[[1]][["res"]]) - ALL_RES[[1]][["res"]][[1]]
for (a in 2:500) {
  P_dot = abind(P_dot,(Reduce(`+`, ALL_RES[[a]][["res"]]) - ALL_RES[[a]][["res"]][[1]]),along = 3)
}
round(Reduce(`+`, ALL_RES[[1]][["res"]]) - ALL_RES[[1]][["res"]][[1]],2)
round(colSums(Reduce(`+`, ALL_RES[[1]][["res"]]) - ALL_RES[[1]][["res"]][[1]]),2)
round(apply(apply(P_dot,c(2,3),sum),c(1),quantile,0.025),3)
round(apply(P_dot,c(1,2),quantile,0.025),3)
round(apply(P_dot,c(1,2),quantile,0.975),3)

round(apply(l_dot,c(1,2),mean)+apply(P_dot,c(1,2),mean),3)


round(apply(l_dot+P_dot-e_dot,c(1,2),mean),5)
round(apply(l_dot+P_dot-e_dot,c(1,2),quantile,0.025),5)
round(apply(l_dot+P_dot-e_dot,c(1,2),quantile,0.5),10)
round(apply(l_dot+P_dot-e_dot,c(1,2),quantile,0.975),5)
(apply(l_dot+P_dot-e_dot,c(1,2),quantile,0.025)*apply(l_dot+P_dot-e_dot,c(1,2),quantile,0.975))<=0 # TRUE is good


p_dot = Reduce(`+`,ALL_RES[[1]][["d_p"]])
for (a in 2:500) {
  p_dot = abind(p_dot,Reduce(`+`,ALL_RES[[a]][["d_p"]]) ,along = 4)
}
round(Reduce(`+`,ALL_RES[[1]][["d_p"]]),2)
round(apply(p_dot,c(1,2,3),quantile,0.025),2)
round(apply(p_dot,c(1,2,3),quantile,0.975),2)
apply(p_dot,c(1,2,3),quantile,0.025) * apply(p_dot,c(1,2,3),quantile,0.975) > 0
# round(apply(Reduce(`+`,ALL_RES[[1]][["d_p"]]), c(1,2), sum),3)
# round(Reduce(`+`, ALL_RES[[1]][["res"]]) - ALL_RES[[1]][["res"]][[1]],3)

p_dot_sum = apply(Reduce(`+`,ALL_RES[[1]][["d_p"]]), c(2,3), sum)
for (a in 2:500) {
  p_dot_sum = abind(p_dot_sum,apply(Reduce(`+`,ALL_RES[[a]][["d_p"]]), c(2,3), sum) ,along = 3)
}
round(apply(Reduce(`+`,ALL_RES[[1]][["d_p"]]), c(2,3), sum),2)
round(apply(p_dot_sum,c(1,2),quantile,0.025),3)
round(apply(p_dot_sum,c(1,2),quantile,0.975),3)
apply(p_dot_sum,c(1,2),quantile,0.025) * apply(p_dot_sum,c(1,2),quantile,0.975) > 0

