library(doFuture)
library(doRNG)
library(coda)
library(rubias)


rubias_MultChains <- function(reference = ref,
                             mixture = mix,
                             nchains = ref %>% select(repunit) %>% unique() %>% unlist() %>% length(),
                             MCMCsteps = 100000,
                             BurnIn = MCMCreps/2
){

MCMCreps <<- MCMCsteps  #saving things as a global variable is typically a no-no :( 
  
#Make variables needed below
  Chn2Samp <- sample(x=seq(1,nchains,by=1),
                     1,
                     replace = F) #which chain do we want to pull results from?
  
### Starting Proportions ###
Baseline_Reps<-reference %>%
    group_by(repunit) %>%
    select(repunit,collection) %>%
    unique()
  
pi_init_List<-lapply(1:length(unique(Baseline_Reps$repunit)),function(x)
    Baseline_Reps %>%
      mutate(RepPi = ifelse(repunit==unique(Baseline_Reps$repunit)[x],1,2))%>%
      mutate(pi_init = ifelse(RepPi==1,
                              0.95/(Baseline_Reps %>%
                                      filter(repunit==unique(Baseline_Reps$repunit)[x]) %>%
                                      nrow()),
                              0.05/(Baseline_Reps %>%
                                      filter(repunit!=unique(Baseline_Reps$repunit)[x]) %>%
                                      nrow())) ) %>%
      ungroup() %>%
      select(collection,pi_init))
  
#Subset on the number of chains, if you want to do fewer than the number of 
#reporting groups just choose nChains randomly
Bl_Chains<-seq(1,length(pi_init_List),by=1) %>% 
  sample(.,size=nchains,replace = F)

pi_init_List<-pi_init_List[Bl_Chains]

### Prior on stock proportions GCg ### 
  # A 1/k prior on baseline populations puts a larger probability on 
  # reporting groups with many populations in the baseline. So we use a 
  # 1/GCg prior to have a more uniform prior over the reporting groups 

#We could give the user the option to change this
  
RepColls<-unique(reference[,c(2,3)])
  
  G<-length(unique(RepColls$repunit))
  PriorRep<-RepColls%>%
    group_by(repunit)%>%
    summarise(nPops=length(unique(collection)))%>%
    mutate(GCg=nPops*G)%>%
    mutate(pi=1/GCg)
  
  prior_GCg<-RepColls%>%
    mutate(pi_param=as.numeric(plyr::mapvalues(x=repunit,
                                               from=PriorRep$repunit,
                                               to=PriorRep$pi,)))%>%
    select(collection,pi_param)

### Run nChains of rubias ### 
registerDoFuture()
plan(multisession) #multisession on windows

X <- 1:nchains

mix_est_list <- foreach(x = X) %dorng% { #dorng vs dopar
  mix_est <- infer_mixture(reference = reference, 
                           mixture = mixture, 
                           gen_start_col = 5,
                           method="MCMC",
                           reps = MCMCreps, 
                           burn_in = BurnIn,
                           pi_init = pi_init_List[[x]],                          
                           pi_prior = prior_GCg)
  mix_est
  
}

plan(sequential)

MCMC_chains<-lapply(1:nchains,function(x)
  mix_est_list[[x]]$mix_prop_traces %>%
    filter(sweep > BurnIn) %>%
    group_by(sweep, repunit) %>%
    summarise(repprop = sum(pi))%>%
    reshape2::dcast(.,formula = sweep ~ repunit , value.var = "repprop") %>%
    select(-sweep)%>%
    coda::mcmc(.)
)%>%
  coda::mcmc.list()

GR_diag <- coda::gelman.diag(MCMC_chains,
                             autoburnin = F, 
                             multivariate = F)

GR_Mat_List<-GR_diag$psrf %>%
  as.data.frame() %>%
  mutate(repunit = row.names(.))

RG_mix_ests <- mix_est_list[[Chn2Samp]]$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) 

trace_subset <- mix_est_list[[Chn2Samp]]$mix_prop_traces %>%
  filter(sweep > BurnIn) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) 

CI_rub <- trace_subset %>%
  group_by(repunit) %>%
  summarise(loCI = quantile(repprop, probs = 0.025),
            hiCI = quantile(repprop, probs = 0.975))

RubiasRes<-cbind(RG_mix_ests[,2:3],CI_rub[,2:3])%>%
  left_join(.,GR_Mat_List,by=c("repunit"))

return(RubiasRes)
} #end of function


## equalize-sample-sizes-fct
# make it a function
# leave-one-out
test_loo_ss <- function(test_2col_df){
  
  test_loo_output <- rubias::assess_reference_loo(test_2col_df, gen_start_col = 5, return_indiv_posteriors = T)
  
  # top p of z for each simulated individual
  top_test_loo <- test_loo_output$indiv_posteriors %>%
    group_by(repunit_scenario, collection_scenario, iter, indiv, simulated_repunit, simulated_collection) %>%
    slice_max(., order_by = PofZ) 
  
  tmp <- top_test_loo %>%
    filter(PofZ > 0.90) %>%
    group_by(simulated_repunit) %>%
    tally(name = "total_per_repunit")
  
  plot_df <- top_test_loo %>%
    filter(PofZ > 0.90) %>%
    group_by(simulated_repunit, repunit) %>%
    tally(name = "assigned_n") %>%
    mutate(correct = ifelse(simulated_repunit == repunit, "TRUE", "FALSE")) %>%
    left_join(., tmp, by = "simulated_repunit") %>%
    mutate(perc_correct = assigned_n/total_per_repunit)
  
  # plot that up
  plot <- plot_df %>%
    ggplot(aes(y = repunit, x = simulated_repunit, size = perc_correct, color = correct)) +
    geom_point() +
    theme_bw() +
    labs(size = "Proportion",
         color = "Accurate assignment",
         y = "Assigned source population",
         x = "Simulated source population")
  
  print(plot)
  
}

### Run Test of function with toy dataset
 # TestOut<-rubias_MultChains(reference = small_chinook_ref,
 #                               mixture = small_chinook_mix,
 #                               nchains = 6,
 #                               MCMCsteps = 1000,
 #                               BurnIn = MCMCreps/2)
