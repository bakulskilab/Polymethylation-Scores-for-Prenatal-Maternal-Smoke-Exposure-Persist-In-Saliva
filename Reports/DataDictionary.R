FFmeta<-read.csv("/scratch/bfoxman_root/bfoxman/blostein/FragileFamilies/Files/FFMetadata_v07.csv")
idvars<-c("id", "ff_id")
demovars<-c("cm1age", "cm1bsex", "cm1ethrace", "m1h3", "m1h3a", "m1h3b", "m1i4", "m1i4a", "m1i4b", "cf1age", "cf1ethrace", "f1h3", "f1h3a", "f1h3b")
smokingvars<-c("m1g4", "f1g4", "f2j5", "f2j5a", "f2j7", "f2j7a", "f3j31", "f3j32", "f4j18", "f4j19", "f5g17", "f5g18", "k5f1k", "k5f1l", "k6d40", "k6d41", "k6d41", "k6d42", "k6d43", "k6d45", "k6d46", "k6d47", "m2j5", "m2j5a", "m2j7", "m2j7a", "m4j18", "m4j19", "m5g17", "m5g18", "n5f17", "n5f18", "p3a22", "p3a23", "p3a23a", "p3a23b", "p3a24", "p4a22", "p4a23", "p5h15", "p5h15b", "p5q3cr", "p6h74", "p6h75", "p6h76", "p6h77", "p6h78")
prenataldrugusevar<-c("m1g2", "m1g3", "m1g5", "m1g6")
FFmyMeta<-FFmeta %>% filter(new_name %in% c(idvars, demovars, smokingvars, prenataldrugusevar))

table(FFmyMeta$type)
myContinuousvars<-FFmyMeta %>% filter(type=="Continuous") %>%pull(new_name)
myBinary<-FFmyMeta%>%filter(type=="Binary")%>%pull(new_name)
myCategorical_O<-FFmyMeta%>%filter(type=="Ordered Categorical")%>%pull(new_name)
myCategorical_U<-FFmyMeta%>%filter(type=="Unordered Categorical")%>%pull(new_name)

