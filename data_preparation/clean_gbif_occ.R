# Clean occurrences GLC23
require(data.table)
require(dplyr)
require(LaF)
require(ggplot2)
require(maptools)
require(sf)
require(rgbif)
require(rgdal)
data(wrld_simpl)
main ="C:/Users/user/pCloud local/boulot/data/geolifeclef/GLC2023/"
#main="/home/cbotella/pcloud/boulot/data/geolifeclef/GLC2023/"

eva="C:/Users/user/pCloud local/boulot/data/EVA/"
dirWrldPoly = 'C:/Users/user/pCloud local/boulot/data/Australian Acacia Worldwide/'

######
# Functions
######

df_as_model = function(filename,sep = ';'){
  names = colnames(fread( filename, sep = sep,nrows= 1 , header = T ))
  names=gsub(" ",".",names)
  names=gsub("'",".",names)
  expr = 'data.frame('
  for (i in 1:length(names)){
    if (i<length(names)){
      expr = paste(expr,names[i],'=NA,',sep="")
    }else{
      expr = paste(expr,names[i],'=NA)',sep="")
    }
  }
  return( eval(parse(text= expr)) )
}

colStat=function(lafCol){
  begin(lafCol)
  gOn=T
  i=1
  while(gOn){
    print(i)
    extr=next_block(lafCol,nrows=1000000)
    if(length(extr)>0){
      if(i==1){
        tab=data.frame(col=extr,count=1)
      }else{
        tab=rbind(tab,
                  data.frame(col=extr,count=1))    
      }
      tab=tab%>%
        group_by(col)%>%
        summarize(count=sum(count))%>%
        as.data.frame()
      i=i+1
    }else{gOn=F}
  }
  return(tab)
}

datasets = data.frame(
  datasetName=c("ArtPortalen",
         "Waarnemingen.be",
         "Estonian Naturalist Soc.",
         "Masaryk Univ Herbarium",
         "Invazivke Slovenia",
         "iRecord verified",
         'Pl@ntNet observations',
         'Swiss Nat Databank',
         'Observation.org',
         'Inventaire forestier IGN',
         'iNaturalist RG',
         'NOR Species Observation',
         'Pl@ntNet automatic',
         'DEN Environmental Portal',
         'Nat. plant monitoring UK'),
  key=c('38b4c89f-584c-41bb-bd8f-cd1def33e92f',
        '7f5e4129-0717-428e-876a-464fbd5d9a47',
        "f1c4df18-12d6-40cb-ab51-5bb0d7f08d6e",
        "54f946aa-2ca9-4a51-9ee5-011219e0381e",
        "ebf3c079-f88e-4b85-bcc5-f568229e68f3",
        "0a013f89-5381-4578-9d82-5f28fd5f1ef6",
        "7a3679ef-5582-4aaa-81f0-8c2545cafc81",
        "83fdfd3d-3a25-4705-9fbe-3db1d1892b13",
        '8a863029-f435-446a-821e-275f4f641165',
        'e5f16d86-e225-4822-97be-a64ce17079c7',
        '50c9509d-22c7-4a22-a47d-8c48425ef4a7',
        'b124e1e0-4755-430f-9eab-894f25a9b59c',
        '14d5676a-2c54-4f94-9023-1e8dcd822aa0',
        '67fabcac-a638-40a6-9bea-aeca8aced9f1',
        '9a6bdcc9-e017-44ea-9cf9-ff6a87fdb8c2'))

colNames=c('gbifID','identifier',
           'modified','publisher',
           'references','datasetName',
           'occurrenceID','recordedBy',
           'eventDate','year','month',
           'day','decimalLatitude','decimalLongitude',
           'coordinateUncertaintyInMeters','taxonRank',
           'datasetKey','species','geoId')

multiplot <- function(plots=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  numPlots = length(plots)
  print(numPlots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

lonlat.to.3035=function(lonLatMatrix){
  # We use the Lambert azimuthal equal area projection
  # https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
  # centered in Europe (EPSG:3035)
  # https://epsg.io/3035-1571
  # the European Environment Agency recommends its usage for European mapping for statistical analysis and display
  # its preserves areas, and distances not much transformed in Europe
  
  proj_ini="+proj=longlat +datum=WGS84"
  projTarget="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  pts = SpatialPointsDataFrame( coords=lonLatMatrix,
                                data=data.frame(id = 1:dim(lonLatMatrix)[1]),
                                proj4string = CRS(proj_ini))
  pts=spTransform(pts,CRS(projTarget))@coords
  colnames(pts)=c('x_EPSG3035','y_EPSG3035')
  return(pts)
}

######
# 1) Filter higher rank and columns
######

keptCols = c("gbifID",
             "identifier" ,
             "modified",
             "publisher",
             'references',
             'collectionCode',
             "datasetName",
             "basisOfRecord",
             'occurrenceID',
             "recordedBy",
             "previousIdentifications",
             'eventDate',
             'year','month','day',
             'decimalLatitude',
             'decimalLongitude',
             "coordinateUncertaintyInMeters",
             "identificationReferences",
             'identificationVerificationStatus',
             "genus",
             'taxonRank',
             'datasetKey',
             'species')

setwd(main)
file= c('occurrence.txt')
saveName ="occ_colFiltered.csv"
taille_serie = 1000000
tmp=df_as_model(file,sep="\t")
selec.ind = which(colnames(tmp)%in%keptCols)
depasse_pas = T
toSkip = 0
erase = T

while (depasse_pas){
  L = data.frame( fread(file, 
                        sep="\t", 
                        nrows=taille_serie,
                        header=FALSE,
                        skip= 1 + toSkip,
                        quote="",
                        col.names = keptCols,
                        select=selec.ind,
                        na.strings = c("NA")) )
  if (dim(L)[1]<(taille_serie-1)){ depasse_pas=FALSE }
  # occurrence is identified at the species level
  L = L[L$taxonRank%in%c('SPECIES','SUBSPECIES'),,drop=F]
  setwd(main)
  if(erase){
    write.table(L,saveName,sep=";",row.names=F,col.names=T)
    erase=F
  }else{
    write.table(L,saveName,sep=";",append = TRUE,row.names=F,col.names=F)}
  ### Number of rows done
  toSkip = toSkip + taille_serie
  print(paste('Number done:',toSkip))
  gc(reset=T)
}

######
# 2) Filter unvalidated identifications 
# from Waarnemingen.be & Artportalen
# and store unique locations
######

Unval = data.frame(
  name=c("artp",
  "waar",
  "esto"),
  key=c('38b4c89f-584c-41bb-bd8f-cd1def33e92f',
               '7f5e4129-0717-428e-876a-464fbd5d9a47',
               "f1c4df18-12d6-40cb-ab51-5bb0d7f08d6e"))
Unval$key=as.character(Unval$key)
# Artportalen and Waarnemingen.be
# Estonian Naturalist Society 

setwd(main)
file= c('occ_colFiltered.csv')
saveName ="occ2_validated.csv"
geoSaveName="geo_tab.csv"
taille_serie = 1000000
depasse_pas = T
toSkip = 0
savedRows=0
savedPos=0
erase = T
geoTab=data.frame(lon=NA,lat=NA,hash=NA,geoId=NA)
maxId = 0
geoTab=geoTab[-1,,drop=F]
while (depasse_pas){
  L = data.frame( fread(file, 
                        sep=";", 
                        nrows=taille_serie,
                        header=FALSE,
                        skip= 1+toSkip,
                        col.names = keptCols,
                        select=1:length(keptCols),
                        na.strings = c("NA")) )
  if (dim(L)[1]<(taille_serie-1)){ depasse_pas=FALSE }
  # occurrence is identified at the species level
  cd = L$datasetKey%in%Unval$key[Unval$name%in%c('waar','artp')]
  validated= cd &
    !is.na(L$identificationVerificationStatus) &
    !L$identificationVerificationStatus%in%c(
      'Not able to validate',
      'Unvalidated',
      'Documentation requested',
      'Dialogue with reporter',
      'Dialogue with validator',
      'unverified')
  Lsave = L[!cd | validated,
            !colnames(L)%in%c(
              'identificationVerificationStatus',
              'previousIdentifications',
              'identificationReferences',
              'genus',
              "basisOfRecord",
              "collectionCode")]
  # Geo id 
  Lsave$hash = paste(Lsave$decimalLongitude,Lsave$decimalLatitude)
  geoTmp=unique(Lsave[,c('decimalLongitude','decimalLatitude','hash')])
  colnames(geoTmp)=c("lon",'lat',"hash")
  geoTmp=geoTmp[!geoTmp$hash%in%geoTab$hash,]
  
  # Saves
  setwd(main)
  if(dim(geoTmp)[1]>0){
    geoTmp$geoId=(maxId+1):(maxId+dim(geoTmp)[1])
    geoTab=rbind(geoTab,geoTmp)
    write.table(geoTab,geoSaveName,sep=";",row.names=F,col.names=T)
  }
  Lsave = merge(Lsave,geoTab[,c('hash','geoId')],by="hash",all.x=T)
  Lsave=Lsave[,!colnames(Lsave)%in%c('hash')]
  
  if(erase){
    write.table(Lsave,saveName,sep=";",row.names=F,col.names=T)
    erase=F
  }else{
    write.table(Lsave,saveName,sep=";",append = TRUE,row.names=F,col.names=F)}
  ### Number of rows done
  toSkip = toSkip + taille_serie
  savedRows=savedRows+dim(Lsave)[1]
  savedPos=dim(geoTab)[1]
  print(paste('Number done:',toSkip))
  print(paste('Total occ number:',savedRows))
  print(paste('Total pos number:',savedPos))
  gc(reset=T)
}

### Filter years prior 2017 and Estonian dataset
file='occ2_validated.csv'
laf=laf_open_csv(file,
                 column_types=c('double','string','string',
                                'string','string','string',
                                'string','string','string',
                                'double','string','string',
                                'double','double','double',
                                'string','string','string','double'),
                 column_names = colNames,sep = ";",
                 dec = ".",skip = 1)

minYear=2017
maxYear=2021
estoKey=datasets$key[datasets$datasetName=="Estonian Naturalist Soc."]
PAkey = datasets$key[datasets$datasetName%in%c("Inventaire forestier IGN",'Nat. plant monitoring UK')]
chunk = 1000000
begin(laf)
gOn=T
i=1
while(gOn){
  print(i)
  extr=next_block(laf,nrows=chunk)
  extr=extr[,!colnames(extr)%in%c('datasetName','modified','identifier')]
  extr$year=round(extr$year)
  extr=extr[(extr$year>=minYear) & extr$year<=maxYear & extr$datasetKey!=estoKey,]
  if(dim(extr)[1]>0){
    cdPA=extr$datasetKey%in%PAkey
    if(i==1){
      tab=extr[!cdPA,]
      PA=extr[cdPA,]
    }else{
      tab=rbind(tab,extr[!cdPA,])
      PA=rbind(PA,extr[cdPA,])
    }
    i=i+1
  }else{gOn=F}
}

setwd(main)
write.table(tab,'occ3_2017_2021_noPA.csv',sep=";",row.names=F,col.names=T)
write.table(PA,'occ3_sup2016_PA.csv',sep=";",row.names=F,col.names=T)
rm(tab,PA,extr,laf);gc(reset=T)

######
# 3) Filter duplicates
######

setwd(main)
gc(reset=T)
file='occ3_2017_2021_noPA.csv'

laf=laf_open_csv(
  file,
  column_types=c('double',
                 'string','string',
                 'string','string','string',
                 'double','string','string',
                 'double','double','double',
                 'string','string','string','double'),
  column_names = colNames[!colNames%in%c('identifier','modified','datasetName')],
  sep = ";",
  dec = ".",
  skip = 1)

chunk = 500000
begin(laf)
gOn=T
i=1
while(gOn){
  print(i)
  extr=next_block(laf,nrows=chunk)
  if(dim(extr)[1]>0){
    if(i==1){tab=extr}else{tab=rbind(tab,extr)}
    i=i+1
  }else{gOn=F}
  gc(reset=T)
}

tab$dateDay = paste(tab$year,tab$month,tab$day,sep="-")
tab$dupId = paste(tab$decimalLongitude,tab$decimalLatitude,tab$species,tab$dateDay,tab$datasetKey,tab$recordedBy)

tmp2=tab%>%
  group_by(decimalLongitude,decimalLatitude,species,dateDay,datasetKey,recordedBy)%>%
  summarize(count=n())%>%
  as.data.frame()

dup = tmp2[tmp2$count>1,]
dup$dupId = paste(dup$decimalLongitude,
                  dup$decimalLatitude,
                  dup$species,
                  dup$dateDay,
                  dup$datasetKey,
                  dup$recordedBy)

# Number removed
sum(dup$count)
ids=1:dim(dup)[1]
toRemove=NULL
for(i in which(dup$count>2)){
  xx=dup$decimalLongitude[i]
  yy=dup$decimalLatitude[i]
  sp=dup$species[i]
  dd=dup$dateDay[i]
  ds=dup$datasetKey[i]
  ob=dup$recordedBy[i]
  cd= tab$dupId==dup$dupId[i]
  bundle = which(cd)
  toRemove = setdiff(toRemove,bundle)
  toRemove=c(toRemove,sample(bundle,length(bundle)-1))
  if(i/100==round(i/100)){
    cat('\r Processed',round(1000*i/length(ids))/10,'%')
    setwd(main)
    save(toRemove,i,dup,file = 'TMP_save.Rdata')
  }
}

toSave=tab[setdiff(1:dim(tab)[1],toRemove),]
toSave=toSave[,!colnames(toSave)%in%c('dateDay','dupId')]

write.table(toSave,'occ3bis_duplic_filtered_noPA.csv',sep=";",row.names=F,col.names=T)

if(F){
  tmp=occ%>%
    group_by(decimalLongitude,decimalLatitude,species)%>%
    summarize(count=n())%>%
    as.data.frame()
  ids=which(tmp$count>=20)
  rownames(datasets)=datasets$key
  for(id in ids){
    sp=tmp$species[id]
    xx=tmp$decimalLongitude[id]
    yy=tmp$decimalLatitude[id]
    cd=occ$species==sp & occ$decimalLongitude==xx & occ$decimalLatitude==yy
    dss=datasets[unique(occ$datasetKey[cd]),'name']
    days=unique(paste(occ$year[cd],occ$month[cd],occ$day[cd]))
    months=unique(occ$month[cd])
    cat('\n ',tmp$count[id],' duplicates for ',sp,' at (',xx,',',yy,') \n ')
    cat(' \n n datasets:',length(dss),'; n obs:',  length(unique(occ$recordedBy[cd])),'; n NA obs:',sum(is.na(occ$recordedBy[cd])),'\n ')
    cat(' \n ',paste(dss,collapse = ","),' \n ')
    cat(' \n n days:',length(days),' \n')
    cat(' \n n months:',length(months),' \n')
  }
}

#######
# 4) count occurrences per species 
#######

gc(reset=T)
setwd(main)
file='occ3bis_duplic_filtered_noPA.csv'
tmp=df_as_model(file)
print(data.frame(colnames(tmp)))
laf=laf_open_csv(
  file,
  column_types=c('double',
                 'string','string',
                 'string','string','string',
                 'double','double','double',
                 'double','double','double',
                 'string','string','string','double'),
  column_names = colNames[!colNames%in%c('identifier','modified','datasetName')],
  sep = ";",
  dec = ".",
  skip = 1)
lafSp = laf$species
tmp=read_lines(lafSp,1:15e6)
tmp=data.frame(species=tmp)
tab=tmp%>%
  cbind(count=1)%>%
  group_by(species) %>%
  summarise(count=sum(count))%>%
  arrange(desc(count)) %>%
  as.data.frame()
tab$rank=1:dim(tab)[1]

res=data.frame(maxC=seq(3e3,6.5e3,50),nOcc=NA)
for(maxC in res$maxC){
  print(maxC)
  fil = tab %>%
    mutate(count2=sapply(count,function(cc)min(cc,maxC))) %>%
    as.data.frame()
  res$nOcc[res$maxC==maxC]= sum(fil$count2)
}
maxC=res$maxC[which.min(abs(5e6-res$nOcc))]
fil = tab %>%
  mutate(count2=sapply(count,function(cc)min(cc,maxC))) %>%
  as.data.frame()

setwd(main)
write.table(fil,
  'occ3_summary_OccPerSpecies_noPA.csv',
  sep=";",row.names=F,col.names=T)

### Plot
thres=c(50,10,5,2)
xs=sapply(thres,function(t)max(fil$rank[fil$count2>=t]))
nOccs=data.frame(x=expand.grid(1:2,xs)[,2],
                 y=rep(c(0,maxC),4),
                 minOcc=expand.grid(1:2,as.character(thres))[,2])
nOccs$minOcc=factor(nOccs$minOcc,levels=as.character(thres))

p=ggplot()+
  geom_line(data=fil,aes(x=rank,y=count2))+
  geom_line(data=nOccs,
    aes(x=x,y=y,group=minOcc,colour=minOcc))+
  scale_colour_manual(values=c('chartreuse4','goldenrod','orange3','orangered3'))+
  ylab('Number of occurrences')+
  xlab('Species rank')+
  scale_x_continuous(breaks=seq(0,max(fil$rank),500))+
  scale_y_continuous(breaks=seq(0,maxC,500))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45),
        text=element_text(size=20))

setwd(main)
png('nOcc_per_species_noPA.png',height=700,width=800)
print(p)
dev.off()

#### En ajoutant les PA
if(T){
  setwd(main)
  PA=read.csv('occ3_sup2016_PA.csv',sep=";",header=T)
  fil2=PA%>%
    mutate(count2=1)%>%
    select(species,count2)%>%
    rbind(fil[,c('species','count2')])%>%
    group_by(species)%>%
    summarise(count2=sum(count2))%>%
    arrange(desc(count2))%>%
    as.data.frame()
  fil2$rank=1:dim(fil2)[1]
  
  setwd(main)
  write.table(fil2,
              'occ3_summary_OccPerSpecies_withPA.csv',
              sep=";",row.names=F,col.names=T)
  
  thres=c(50,10,5,2)
  xs=sapply(thres,function(t)max(fil2$rank[fil2$count2>=t]))
  nOccs=data.frame(x=expand.grid(1:2,xs)[,2],
                   y=rep(c(0,maxC),4),
                   minOcc=expand.grid(1:2,as.character(thres))[,2])
  nOccs$minOcc=factor(nOccs$minOcc,levels=as.character(thres))
  
  p=ggplot()+
    geom_line(data=fil2,aes(x=rank,y=count2))+
    geom_line(data=nOccs,
              aes(x=x,y=y,group=minOcc,colour=minOcc))+
    scale_colour_manual(values=c('chartreuse4','goldenrod','orange3','orangered3'))+
    ylab('Number of occurrences')+
    xlab('Species rank')+
    scale_x_continuous(breaks=seq(0,max(fil2$rank),500))+
    scale_y_continuous(breaks=seq(0,max(fil2$count2),500))+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45),
          text=element_text(size=20))
  
  setwd(main)
  png('nOcc_per_species_withPA.png',height=700,width=800)
  print(p)
  dev.off()
  
  thres=c(50,10,5,2)
  xs=sapply(thres,function(t)max(fil2$rank[fil2$count2>=t]))
  nOccs=data.frame(x=expand.grid(1:2,xs)[,2],
                   y=rep(c(0,log10(max(fil2$count2))),4),
                   minOcc=expand.grid(1:2,as.character(thres))[,2])
  nOccs$minOcc=factor(nOccs$minOcc,levels=as.character(thres))
  
  plog=ggplot()+
    geom_line(data=fil2,aes(x=rank,y=log10(count2)))+
    geom_line(data=nOccs,
              aes(x=x,y=y,group=minOcc,colour=minOcc))+
    scale_colour_manual(values=c('chartreuse4','goldenrod','orange3','orangered3'))+
    ylab('Number of occurrences')+
    xlab('Species rank')+
    scale_x_continuous(breaks=seq(0,max(fil2$rank),500))+
    scale_y_continuous(breaks=seq(0,round(max(log10(fil2$count2)))+1,0.5))+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45),
          text=element_text(size=20))
  
  setwd(main)
  png('nOcc_per_species_withPA_log10.png',height=700,width=800)
  print(plog)
  dev.off()
}


#######
# 5) Filter & Subsample to 5M occ
#######

gc(reset=T)
file='occ3bis_duplic_filtered_noPA.csv'
setwd(main)
sps=read.csv('occ3_summary_OccPerSpecies_noPA.csv',sep=";",header=T)
colNames=colnames(df_as_model(file))
laf=laf_open_csv(
  file,
  column_types=c('double',
                 'string','string',
                 'string','string','string',
                 'double','string','string',
                 'double','double','double',
                 'string','string','string','double'),
  column_names = colNames,
  sep = ";",
  dec = ".",
  skip = 1)

chunk = 1000000
begin(laf)
gOn=T
i=1
while(gOn){
  print(i)
  extr=next_block(laf,nrows=chunk)
  if(dim(extr)[1]>0){
    if(i==1){tab=extr}else{tab=rbind(tab,extr)}
    i=i+1
  }else{gOn=F}
  gc(reset=T)
}

maxC=max(sps$count2)
# List of subsampled species 
spVec=sps$species[sps$count2==maxC]
# Initialize empty table where to put subsampled occ 
finalTab = tab[1:sum(sps$count2),]
finalTab[]=NA
# Begin table with non-subsampled species
cd=!tab$species%in%spVec
starto=sum(cd)
finalTab[1:sum(cd),]=tab[cd,]
# For each species to subsample
for(sp in spVec){
  tmp = tab[tab$species==sp,]
  # Subsample randomly maxC occurrences
  tmp = tmp[sample(1:dim(tmp)[1],maxC),]
  # Add to finalTab
  finalTab[(starto+1):(starto+maxC),]=tmp
  starto=starto+maxC
  gc(reset=T)
  cat('\r Processed ',round(1000*which(spVec==sp)/length(spVec))/10,'%')
}
finalTab=finalTab[!is.na(finalTab$gbifID),]
# Randomly SHUFFLE finalTab
finalTab=finalTab[sample(1:dim(finalTab)[1],dim(finalTab)[1]),]

setwd(main)
write.table(finalTab,'occ4_filtMaxOcc.csv',sep=";",col.names=T,row.names=F)

table(finalTab$datasetKey)%>%
  as.numeric()%>%
  cbind(key=names(table(finalTab$datasetKey)))%>%
  merge(datasets,by="key",all.x=T)

#######
# 6) Make shared PA/PO CSV files 
#######

### Prepare PO file
setwd(main)
occ=read.csv('occ4_filtMaxOcc.csv',sep=";",header=T)
colnames(occ)[c(10,11,12,5)]=c('lat','lon','geoUncertaintyInM','observer')
occ=occ%>%
  mutate(date=substr(eventDate,1,10),
         dayOfYear=as.numeric(strftime(eventDate,format='%j')))%>%
  merge(datasets,by.x='datasetKey',by.y="key",all.x=T)

#select(lon,lat,x_EPSG3035,y_EPSG3035,geoUncertaintyInM,
#         year,dayOfYear,date,datasetName,observer,species)

maxID=dim(occ)[1]
# Shuffle
occ=occ[sample(1:maxID,maxID),]

PO=cbind(occ,lonlat.to.3035(occ[,c("lon","lat")]))%>%
  mutate(glcID=1:maxID)%>%
  select(glcID,gbifID,lon,lat,x_EPSG3035,y_EPSG3035,
       geoUncertaintyInM,
       year,dayOfYear,date,
       datasetName,observer,
       species)
#PO=read.csv('PO_gbif.csv',sep=";",header=T)
PO=PO[PO$lon>(-13) | PO$lat>=40.5,]

write.table(PO,'PO_gbif.csv',sep=";",row.names=F,col.names=T)


### Prepare PA_gbif file
setwd(main)
PA=read.csv('occ3_sup2016_PA.csv',sep=";",header=T)
colnames(PA)[c(10:12,5)]=c('lat','lon',"geoUncertaintyInM","observer")

# Get dataset name
PA=PA%>%
  merge(datasets,by.x="datasetKey",by.y="key",all.x=T)
PA$dayOfYear=as.numeric(strftime(PA$eventDate,format='%j'))

# Aggregate releves per x,y,year,datasetKey
sites=PA%>%
  group_by(lon,lat,year,datasetName)%>%
  summarise(nDates=length(unique(eventDate)),
            richness=length(unique((species))))%>%
  as.data.frame()
# Remove plots that had only one releve in a year 
#sites=sites[sites$datasetName!='Nat. plant monitoring UK' | sites$nDates>1,]

# Set PlotID for GBIF plots
sites=sites%>%
  group_by(lon,lat,year)%>%
  summarise(n=n())
sites$PlotID=3000000+(1:dim(sites)[1])

# Associate PlotID
PA=PA%>%
  merge(sites[,c('lon','lat','year','PlotID')],by=c('lon','lat','year'),all.x=T)%>%
  filter(!is.na(PlotID))

nPAgbif=dim(PA)[1]
PA=PA[sample(1:nPAgbif,nPAgbif),]
PA$glcID=(maxID+1):(maxID+nPAgbif)

PAgbif=PA%>%
  group_by(glcID,gbifID,PlotID,lon,lat,year,dayOfYear,datasetName,species,observer)%>%
  summarise(geoUncertaintyInM=min(geoUncertaintyInM,na.rm=T))%>%
  mutate(PlotObservationID_eva=NA,
         date=strptime(paste(substr(year,3,4), dayOfYear), "%y %j"),
         source="GBIF",
         access='public')
PAgbif=cbind(PAgbif,lonlat.to.3035(data.frame(lon=PAgbif$lon,lat=PAgbif$lat)))%>%
  select(glcID,gbifID,PlotID,
         lon,lat,x_EPSG3035,y_EPSG3035,
         geoUncertaintyInM,
         year,dayOfYear,date,
         datasetName,observer,
         species,source,access,PlotObservationID_eva)

setwd(main)
write.table(PAgbif,'PA_gbif.csv',sep=";",row.names=F,col.names=T)

### Prepare EVA datasets 
setwd(eva)
#PAeva=read.csv('occ_PA_eva_CBNA_CBNmed.csv',sep=";",header=T)
stat=read.csv('EVA_datasets_preSelec.csv',sep="\t",header=T)
dsToKeep=stat$Dataset[stat$nPlotSup2017>60]

setwd(eva)
sp=read.csv('153_PlantNet20220701_notJUICE_species.csv',sep="\t",header=T)

tt=read.csv('153_PlantNet20220701_notJUICE_header.csv',sep="\t",header=T)
# Filter and clean dates
tt=tt[nchar(tt$Date.of.recording)==10,]
tt$date = paste0(substr(tt$Date.of.recording,7,10),
                 '-',substr(tt$Date.of.recording,4,5),
                 '-',substr(tt$Date.of.recording,1,2))
tt$year=as.numeric(substr(tt$date,1,4))
tt$dayOfYear=as.numeric(strftime(tt$date,format='%j'))
tt=tt%>%
  select(PlotObservationID,Longitude,
         Latitude,Location.uncertainty..m.,year,
         dayOfYear,Dataset,Author)
colnames(tt)[c(2:4,7:8)]=c('lon','lat','geoUncertaintyInM','datasetName','observer')
tt=tt%>%
  filter(is.na(geoUncertaintyInM) | geoUncertaintyInM<=100 & datasetName%in%dsToKeep)%>%
  group_by(PlotObservationID,lon,lat,year,dayOfYear,datasetName,observer)%>%
  summarise(geoUncertaintyInM=min(geoUncertaintyInM,na.rm=T))
tt$geoUncertaintyInM[is.infinite(tt$geoUncertaintyInM)]=NA
# Set PlotID for EVA plots
plots=tt%>%
  group_by(lon,lat,year)%>%
  summarise()
plots$PlotID=1:dim(plots)[1]
tt=merge(tt,plots,by=c('lon','lat','year'),all.x=T)


# Filter species out of these datasets
cd = sp$PlotObservationID%in%PAev$PlotObservationID & sp$Taxon.group=='Vascular plant' & sp$Cover..>0
spSel=sp[cd,c('PlotObservationID','Taxonomy','Turboveg2.concept','Matched.concept','Match','Original.taxon.concept')]

write.table(spSel,'species_target_datasets.csv',sep=";",col.names=T,row.names=F)
sp=spSel;rm(spSel)
gc(reset=T)
#sp=read.csv('species_CBNa_CBNmed.csv',sep=";",header=T)
sp=sp[sp$Match!=0,]

sps=data.frame(initName=unique(sp$Turboveg2.concept))
sps=sps[!is.na(sps$initName),,drop=F]
# Taxonomic matching: Get GBIF scientificName
cols = c("usageKey","rank","scientificName",'phylum','class','kingdom',
         "matchType","confidence","synonym",
         "status","canonicalName","species","speciesKey")
toAdd=as.data.frame(matrix(NA,dim(sps)[1],length(cols)))
colnames(toAdd)=cols
for(i in 1:dim(sps)[1]){
  tente = as.data.frame(rgbif::name_backbone(as.character(sps$initName[i])))
  if(sum(!cols%in%colnames(tente))>0){
    for(col in cols[!cols%in%colnames(tente)]){
      eval(parse(text=paste('tente$',col,'=NA',sep="")))
    }
  }
  tente = tente[,cols]
  if(!is.null(dim(tente))){
    toAdd[i,]= tente[1,]
  }
  if(round(i/100)==i/100){
    cat('\r Processed...',round(10000*i/dim(sps)[1])/100,'%')
    write.table(toAdd,'eva_species_match_gbif.csv',sep=";",col.names=T,row.names=F)
  }
}
write.table(toAdd,'eva_species_match_gbif.csv',sep=";",col.names=T,row.names=F)
sps_ = cbind(sps,toAdd)

table(sps_$matchType)
table(sps_$rank)
sps_=sps_[sps_$rank%in%c('SPECIES','SUBSPECIES','VARIETY'),]
# Validate fuzzy matches
if(T){
  validated=c('Carex ferruginea subsp. australpina',
              'Ranunculus keupferi',
              'Campanula cochlearifolia',
              'Carex hallerana',
              'Alopecurus gerardi',
              'Erysimum jugicola',
              'Centaurea triumfetti',
              'Arabis soyeri subsp. subcoriacaea',
              'Limonium auriculae-ursifolium',
              'Melilotus indica',
              'Melilotus neapolitana',
              'Euphorbia seguierana subsp. seguierana',
              'Odontites lutea subsp. lutea',
              'Melilotus sulcata',
              'Bidens tripartitus',
              'Euphorbia seguierana',
              'Armeria girardii',
              'Rhamnus catharticus',
              'Centranthus lecoquii',
              'Melica bauhinii',
              'Platycapnos spicatus',
              'Spergularia bocconii',
              'Ranunculus parnassiifolius subsp. parnassiifolius',
              'Ranunculus parnassiifolius subsp. heterocarpus',
              'Ranunculus parnassiifolius',
              'Carex mairii',
              'Centaurea spinabadia subsp. hanryi',
              'Ranunculus gouanii',
              'Aira tenorii',
              'Ampelodesmos mauritanica',
              'Erysimum montosicola',
              'Trochiscanthes nodiflora',
              'Centaurea triumfetti subsp. aligera',
              'Hymenolobus procumbens var. revelieri',
              'Gnaphalium luteo-album',
              'Bromus hordeaceus subsp. thominii',
              'Carduncellus monspelliensium',
              'Centaurea triumfetti subsp. semidecurrens',
              'Odontites viscosa',
              'Ranunculus revelieri',
              'Odontites lutea',
              'Taraxacum braunblanquetii',
              'Galium pseudaristatum',
              'Odontites lanceolata',
              'Juncus acutus subsp. leopoldii')
}
sps_$matchType[sps_$initName%in%validated]="EXACT"
# Remove remaining fuzzy matches
sps_ = sps_[sps_$matchType%in%c("EXACT",'HIGHERRANK'),c('initName','species','speciesKey','confidence')]
colnames(sps_)=c('name_eva','species','speciesKey_GBIF','match_confidence')

spli=sp[sp$Turboveg2.concept%in%sps_$name_eva,]%>%
  merge(sps_,by.x='Turboveg2.concept',by.y='name_eva',all.x=T)
spli=spli[,c('PlotObservationID','Turboveg2.concept','Match','species')]

setwd(eva)
write.table(sps_,'eva_species_list_match_gbif.csv',sep=";",row.names=F,col.names=T)
write.table(spli,'species_vs_plotObs_filtered.csv',sep=";",row.names=F,col.names=T)

spli=read.csv('species_vs_plotObs_filtered.csv',sep=";",header=T)

# Merge with species 
PAev=tt%>%
  merge(spli[,c('PlotObservationID','species')],by="PlotObservationID",all.x=T)%>%
  mutate(date=strptime(paste(substr(year,3,4), dayOfYear), "%y %j"),
         source="EVA",
         access=ifelse(datasetName%in%c("CBNMed","CBNA"),'public','private'),
         PlotObservationID_eva=PlotObservationID)
PAev=cbind(PAev,lonlat.to.3035(data.frame(lon=PAev$lon,lat=PAev$lat)))%>%
  select(PlotID,
         lon,lat,x_EPSG3035,y_EPSG3035,
         geoUncertaintyInM,
         year,dayOfYear,date,
         datasetName,observer,
         species,source,access,PlotObservationID_eva)

setwd(main)
write.table(PAev,'PA_EVA.csv',sep=";",row.names=F,col.names=T)
rm(tt,spli);gc(reset=T)

### Assemble PA from EVA and GBIF
PAev=read.csv('PA_EVA.csv',sep=";",header=T)
PAgbif=read.csv('PA_gbif.csv',sep=";",header=T)

public=PAev[PAev$access=="public",]
public$gbifID=NA
maxID=max(PAgbif$glcID)
nPApu=dim(public)[1]
public=public[sample(1:nPApu,nPApu),]
public=public%>%
  mutate(glcID=(maxID+1):(maxID+nPApu),
         date=as.POSIXct(date))%>%
  select(glcID,gbifID,PlotID,lon,lat,x_EPSG3035,y_EPSG3035,geoUncertaintyInM,year,dayOfYear,date,datasetName,observer,species,source,access,PlotObservationID_eva)

PAdispo=rbind(PAgbif,public)
PAdispo=PAdispo%>%
  filter(year>=2017 & year<=2021)%>%
  filter(species%in%unique(PO$species))

write.table(PAdispo[,colnames(PAdispo)!="access"],'PA_public.csv',sep=";",row.names=F,col.names=T)

PApriv=PAev[PAev$access=="private",colnames(PAev)!="access"]

write.table(PApriv,'PA_private.csv',sep=";",row.names=F,col.names=T)

# J'ai donc ajouté les fichiers mentionnés hier, PO_public (presence-only, pour le train), PA_public (presence-absence, pour train / test, contient IGN, NPMS, CBNMed, CBNA), et enfin la face cachée de l'iceberg EVA PA_private 
# le PlotID correspond aux combinaisons de lon,lat,year
#Attention, il peut y avoir plusieurs fois la espèce pour un PlotID (ex: relevée à différentes dates, observateurs différents)
#Pour sentinel/bioclim, la combinaison PlotID-dayOfyear sera l'unité d'extraction et l'identifiant, faut pas extraire pour chaque ligne
#Pour le split train / test, voir les colonnes x_EPSG3035 y_EPSG3035 

# N.B. taxonomie d'EVA: J'ai dû matcher tous les noms EVA avec le GBIF backbone taxonomy, le nom d'espèce obtenu est la colonne "species". 
# Quand ça matchait pas ou au dessus du rang espèce, j'ai viré. Quand le matching était approximatif, SI c'était CBN j'ai validé à la mano
# mais pour les autres datasets, j'ai tout viré mais c'était des taxons peu fréquents (<1.2% de tous les plots EVA). 
# En tout j'ai dû virer environ 300 noms d'EVA sur >13,000

######
# 6)bis Anonymise species
######

# After split, check no train species absent from test, if so remove from test
# Add readme for each predictor folder with citation 

setwd(main)
PA=read.csv('PA_public.csv',sep=";",header=T)
PO=read.csv('PO_gbif.csv',sep=";",header=T)

#length(setdiff(unique(PA$species),unique(PO$species)))
set.seed(32)
# All GLC species are present in PO 
spNames=data.frame(species=unique(PO$species))
spNames=spNames[sample(1:dim(spNames)[1],dim(spNames)[1]),,drop=F]
spNames$speciesId=1:dim(spNames)[1]

PO_ano=PO%>%
  merge(spNames,by="species",all.x=T)
PO_ano=PO_ano[,colnames(PO_ano)!='species']

# Set PlotID for PO
pIds=PO_ano%>%
  group_by(lon,lat,year)%>%
  summarise(n=n())
pIds$PlotID=(max(PA$PlotID)+1):(max(PA$PlotID)+dim(pIds)[1])
PO_ano=PO_ano%>%
  merge(pIds[,c("lon","lat",'year',"PlotID")],by=c("lon","lat",'year'),all.x=T)

# Filter duplicates between PL obs and queries
obs=PO_ano[PO_ano$datasetName=="Pl@ntNet observations",]%>%
  mutate(checkID=paste0(PlotID,'_',dayOfYear,'_',observer,'_',speciesId))
rest=PO_ano[PO_ano$datasetName!="Pl@ntNet observations",]%>%
  mutate(checkID=paste0(PlotID,'_',dayOfYear,'_',observer,'_',speciesId))
PO_ano=PO_ano[!PO_ano$glcID%in%obs$glcID[obs$checkID%in%rest$checkID],]

# Filter duplicates with same PlotID,dayOfYear,speciesId,observer,datasetName
fir=function(vec){vec[1]}
PO_ano = PO_ano%>%
  group_by(PlotID,dayOfYear,speciesId,observer,datasetName)%>%
  summarise(glcID=fir(glcID),
            gbifID=fir(gbifID),
            date=fir(date),
            lon=fir(lon),
            lat=fir(lat),
            x_EPSG3035=fir(x_EPSG3035),
            y_EPSG3035=fir(y_EPSG3035),
            geoUncertaintyInM=min(geoUncertaintyInM),
            year=fir(year))

PA_ano=PA%>%
  merge(spNames,by="species",all.x=T)
PA_ano=PA_ano[,colnames(PA_ano)!='species']

setwd(main)
write.table(spNames,'spNamesTable.csv',sep=";",row.names=F,col.names=T)
write.table(PO_ano,'PO_anonymised.csv',sep=";",row.names=F,col.names=T)
write.table(PA_ano,"PA_public_anonymised.csv",sep=";",row.names=F,col.names=T)

#######
# 7) PA sites analysis
#######

data(wrld_simpl)
poly = st_as_sf(wrld_simpl)
st_crs(poly) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
borders=poly[poly$NAME%in%c("France","United Kingdom"),]

setwd(main)
PA=read.csv('PA_public.csv',sep=";",header=T)

dir="C:/Users/user/pCloud local/boulot/"

sites=PA%>%
  group_by(lon,lat,year,datasetName)%>%
  summarise(nDates=length(unique(dayOfYear)),
            richness=length(unique((species))))

p=ggplot(sites,aes(x=richness,colour=datasetName))+
  geom_histogram(position="identity",bins = 30,alpha=.01,size=1)+theme_bw()+
  scale_x_continuous(breaks=seq(0,100,5))
print(p)

setwd(main)
png('PA_richness_per_dataset.png',width=700)
print(p)
dev.off()

p2=ggplot(sites[sites$dataset=='NPMS UK',],aes(x=nDates,group=nDates,y=richness))+
  geom_boxplot()+
  geom_violin(alpha=.4)+
  scale_y_continuous(breaks=seq(0,80,5))

setwd(main)
png('PA_richness_vs_nDates.png',width=700)
print(p2)
dev.off()

cd=sites$dataset=='Forest inventory IGN' | sites$nDates>=2
p3=ggplot()+
  geom_smooth(data=sites[cd,],aes(x=lat,y=richness,group=dataset,color=dataset),
              size=1)+
  scale_color_manual(values=c('darkblue','darkred'))+
  geom_point(data=sites[cd,],aes(x=lat,y=richness,colour=dataset),alpha=.3,pch=15)+
  xlab('Latitude')+
  theme_bw()+theme(text=element_text(size=25))

setwd(main)
png('PA_richness_vs_latitude_NPMSwellSamp.png',height=900,width=1100)
print(p3)
dev.off()

###
p4=ggplot()+
  geom_smooth(data=sites,aes(x=lat,y=richness,group=dataset,color=dataset),
              size=1)+
  scale_color_manual(values=c('darkblue','darkred'))+
  geom_point(data=sites[cd,],aes(x=lat,y=richness,colour=dataset),alpha=.3,pch=15)+
  xlab('Latitude')+
  theme_bw()+theme(text=element_text(size=25))

setwd(main)
png('PA_richness_vs_latitude_allSites.png',height=900,width=1100)
print(p4)
dev.off()

#### CARTE
p5=ggplot()+geom_sf(data=borders)+
  geom_point(data=sites,aes(x=lon,y=lat,colour=richness),size=1,pch=15,alpha=.4)+
  scale_colour_gradient2(mid="darkorchid4",high="goldenrod")+
  xlab('Longitude')+ylab('Latitude')+
  theme_bw()+theme(text=element_text(size=20))

setwd(main)
png('PA_richness_map.png',height=2000,width=1600)
print(p5)
dev.off()

#### CARTE PA richesse
p5=ggplot()+geom_sf(data=borders)+
  geom_point(data=sites,aes(x=lon,y=lat,colour=richness),size=.8,pch=15,alpha=.4)+
  scale_colour_gradient2(mid="darkorchid4",high="goldenrod")+
  xlab('Longitude')+ylab('Latitude')+
  theme_bw()+theme(text=element_text(size=20))

setwd(main)
png('PA_richness_map.png',height=2000,width=1600)
print(p5)
dev.off()

#### CARTE PA par dataset
p5=ggplot()+geom_sf(data=borders)+
  geom_point(data=sites,aes(x=lon,y=lat,colour=datasetName),size=.7,pch=15,alpha=.4)+
  xlab('Longitude')+ylab('Latitude')+
  theme_bw()+theme(text=element_text(size=20))

setwd(main)
png('PA_dataset_map.png',height=1500,width=1400)
print(p5)
dev.off()



#### Pertes SI <2019
# Perte de 3000 releves sur 19000
relev=PA%>%
  group_by(lon,lat,year,datasetKey)%>%
  summarise(nDates=length(unique(eventDate)),
            richness=length(unique((species))))%>%
  merge(ds,by="datasetKey",all.x=T)%>%
  as.data.frame()
relev$relId=1:dim(relev)[1]

print(dim(relev)[1])
print(sum(relev$year<=2019))

toKeep = unique(relev$relId[relev$year<=2019])

# n Releves
PA=PA%>%
  merge(relev[,c('lon','lat','year','datasetKey','dataset','relId')],
        by=c('lon','lat','year','datasetKey'),
        all.x=T)

### Perte de 700 sites géographiques sur 14400
sitesAll=PA%>%
  group_by(lon,lat,dataset)%>%
  summarise(nRel=length(unique(relId)),
            nDates=length(unique(eventDate)),
            richness=length(unique((species))))%>%
  as.data.frame()
print(dim(sitesAll)[1])

sitesBe2019=PA%>%
  filter(relId%in%toKeep)%>%
  group_by(lon,lat,dataset)%>%
  summarise(nRel=length(unique(relId)),
            nDates=length(unique(eventDate)),
            richness=length(unique((species))))%>%
  as.data.frame()
print(dim(sitesBe2019)[1])  

### Perte de 80 espèces
length(unique(PA$species))
length(unique(PA$species[PA$relId%in%toKeep]))



######
# Ground Truth for Kaggle
######

require(dplyr)
setwd(main)
train=read.csv('PA_public_anonymised_train.csv',sep=";",header=T)
test=read.csv('PA_public_anonymised_test.csv',sep=";",header=T)

gt=test%>%
  group_by(PlotID,dayOfYear)%>%
  summarise(Expected=paste(sort(unique(speciesId)),collapse=" "))
gt$Id=1:dim(gt)[1]

# Updated PA_public_anonymised_test.csv
tmp=test%>%
  merge(gt[,c('PlotID','dayOfYear','Id')],by=c('PlotID','dayOfYear'),all.x=T)

# Ground truth
gt=gt[,c('Id','Expected')]

# Perfect submission file (no error)
subPerf=gt%>%
  rename(Predicted=Expected)

# Baseline submission (based on train PA)
toCal=train%>%
  group_by(PlotID,dayOfYear,speciesId)%>%
  summarise()
freqPerSp=toCal%>%
  group_by(speciesId)%>%
  summarise(prop=n()/dim(.)[1])

calib=data.frame(minF=seq(1e-4,.03,5e-4),nSp=NA,score=NA,testScore=NA)
for(minF in calib$minF){
  pred=freqPerSp$speciesId[freqPerSp$prop>minF]
  scor=toCal%>%
    group_by(PlotID,dayOfYear)%>%
    summarise(tps=length(base::intersect(speciesId,pred)),
              fps=length(base::setdiff(pred,speciesId)),
              fns=length(base::setdiff(speciesId,pred)))%>%
    mutate(f1=tps/(tps+(fps+fns)/2))
  Ttps=sapply(GT,function(truth)length(intersect(pred,as.numeric(truth))))
  Tfps=sapply(GT,function(truth)length(setdiff(pred,as.numeric(truth))))
  Tfns=sapply(GT,function(truth)length(setdiff(as.numeric(truth),pred)))
  
  calib$nSp[calib$minF==minF]=length(pred)
  calib$testScore[calib$minF==minF]=mean(Ttps/(Ttps+(Tfps+Tfns)/2))
  calib$score[calib$minF==minF]= mean(scor$f1)
  cat('\n Num sp:',length(pred))
}

png('microF1_vs_minPrevalence.png')
plot(calib$nSp,calib$score,col="red",type="l",xlab = "# included species in constant baseline",ylab="micro F1 score")
lines(calib$nSp,calib$testScore,col='green',type="l")
dev.off()

minF=calib$minF[calib$score==max(calib$score)]
pred=sort(freqPerSp$speciesId[freqPerSp$prop>=minF])
print(length(pred))
print(mean(sapply(GT,function(truth)length(as.numeric(truth)))))
print(mean((toCal%>%group_by(PlotID,dayOfYear)%>%count())$n))

print(paste('micro F score Kaggle attendu=',calib$testScore[calib$score==max(calib$score)]))
baseline=data.frame(Id=gt$Id,Predicted=paste(pred,collapse = " "))


write.table(gt,'gt.csv',sep=",",row.names=F,col.names=T,quote = F)
write.table(subPerf,'sub_perfect.csv',sep=",",row.names=F,col.names=T,quote = F)
write.table(baseline,'baseline.csv',sep=",",row.names=F,col.names=T,quote = F)


write.table(tmp,'PA_public_anonymised_test.csv',sep=";",row.names=F,col.names=T)

######
# Make Test CSV without speciesId
###### 

setwd(main)
test=read.csv('Presences_Absences_test_noId.csv',sep=";",header=T)
witId=read.csv('PA_public_anonymised_test.csv',sep=";",header=T)

testF=test%>%
  merge(witId[,c('glcID','Id')],by="glcID",all.x=T)
testF=testF[,!colnames(testF)%in%c(
  'PlotObservationID_eva',
  "year_ecodatacube_quarter",
  "source")]
blind=testF[,colnames(testF)!="speciesId"]

write.table(testF,'Presences_Absences_test.csv',sep=";",row.names=F,col.names=T)
write.table(blind,'Presences_Absences_test_blind.csv',sep=";",row.names=F,col.names=T)


setwd(main)
blind=read.csv('Presences_Absences_test_blind.csv',sep=";",header=T)
head(blind)
blindag=blind%>%
  group_by(Id,datasetName,date,dayOfYear,year,patchID,timeSerieID)%>%
  summarise(lon=mean(lon),
            lat=mean(lat),
            x_EPSG3035=mean(x_EPSG3035),
            y_EPSG3035=mean(y_EPSG3035),
            geoUncertaintyInM=max(geoUncertaintyInM))


write.table(blindag,'test_blind.csv',sep=";",row.names=F,col.names=T)

#######
# EVA datasets
######


setwd(main)
PA=read.csv('occ3_sup2016_PA.csv',sep=";",header=T)
setwd(eva)
sp=read.csv('153_PlantNet20220701_notJUICE_species.csv',sep="\t",header=T)

tt=read.csv('153_PlantNet20220701_notJUICE_header.csv',sep="\t",header=T)

# Quick summary of EVA plots
tt=tt[nchar(tt$Date.of.recording)==10,]
tt$year=sapply(tt$Date.of.recording,function(stri)strsplit(stri,split = "\\.")[[1]][3])
tt$year=as.numeric(tt$year)

mostComon=function(vec){
  tab=table(vec)
  res=names(tab)[tab==max(tab)][1]
  return(res)
}

sel=tt[tt$year>=2014,]
stat=sel%>%
  group_by(Dataset)%>%
  summarise(nPlotSup2014=n(),
            nPlotSup2016=sum(year>=2016),
            medArea_m2=median(Relevé.area..m..,na.rm=T),
            countryMaj=mostComon(Country),
            for_EVA=sum(FOR_EVA=="Y"),
            xMin=min(Longitude),
            xMax=max(Longitude),
            yMin=min(Latitude),
            yMax=max(Latitude))%>%
  filter(nPlotSup2014>100)%>%
  arrange(desc(nPlotSup2014))%>%
  as.data.frame()

print(stat)

setwd(eva)
write.table(stat,'EVA_datasets_preSelec.csv',sep=";",row.names=F,col.names=T)

#PlotObservationID
# Keep only plots from CBNA & CBNmed
dss=c('CBNMed','CBNA')
sites=tt[tt$Dataset%in%dss,]
cd = sp$PlotObservationID%in%sites$PlotObservationID & sp$Taxon.group=='Vascular plant' & spSel$Cover..>0
spSel=sp[cd,c('PlotObservationID','Taxonomy','Turboveg2.concept','Matched.concept','Match','Original.taxon.concept')]

write.table(spSel,'species_CBNa_CBNmed.csv',sep=";",col.names=T,row.names=F)
sp=spSel;rm(spSel)
gc(reset=T)
#sp=read.csv('species_CBNa_CBNmed.csv',sep=";",header=T)
sp=sp[sp$Match!=0,]

sps=data.frame(initName=unique(sp$Turboveg2.concept))
# Get GBIF scientificName
cols = c("usageKey","rank","scientificName",'phylum','class','kingdom',
         "matchType","confidence","synonym",
         "status","canonicalName","species","speciesKey")
toAdd=as.data.frame(matrix(NA,dim(sps)[1],length(cols)))
colnames(toAdd)=cols
for(i in 1:dim(sps)[1]){
  tente = as.data.frame(rgbif::name_backbone(as.character(sps$initName[i])))
  if(sum(!cols%in%colnames(tente))>0){
    for(col in cols[!cols%in%colnames(tente)]){
      eval(parse(text=paste('tente$',col,'=NA',sep="")))
    }
  }
  tente = tente[,cols]
  if(!is.null(dim(tente))){
    toAdd[i,]= tente[1,]
  }
  if(round(i/100)==i/100){
    cat('\r Processed...',round(10000*i/dim(sps)[1])/100,'%')
    write.table(toAdd,'tmp_species_match_gbif.csv',sep=";",col.names=T,row.names=F)
  }
}
sps_ = cbind(sps,toAdd)

table(sps_$matchType)
table(sps_$rank)
sps[sps_$matchType=="FUZZY",]
# Validate fuzzy matches
if(T){
  validated=c('Carex ferruginea subsp. australpina',
              'Ranunculus keupferi',
              'Campanula cochlearifolia',
              'Carex hallerana',
              'Alopecurus gerardi',
              'Erysimum jugicola',
              'Centaurea triumfetti',
              'Arabis soyeri subsp. subcoriacaea',
              'Limonium auriculae-ursifolium',
              'Melilotus indica',
              'Melilotus neapolitana',
              'Euphorbia seguierana subsp. seguierana',
              'Odontites lutea subsp. lutea',
              'Melilotus sulcata',
              'Bidens tripartitus',
              'Euphorbia seguierana',
              'Armeria girardii',
              'Rhamnus catharticus',
              'Centranthus lecoquii',
              'Melica bauhinii',
              'Platycapnos spicatus',
              'Spergularia bocconii',
              'Ranunculus parnassiifolius subsp. parnassiifolius',
              'Ranunculus parnassiifolius subsp. heterocarpus',
              'Ranunculus parnassiifolius',
              'Carex mairii',
              'Centaurea spinabadia subsp. hanryi',
              'Ranunculus gouanii',
              'Aira tenorii',
              'Ampelodesmos mauritanica',
              'Erysimum montosicola',
              'Trochiscanthes nodiflora',
              'Centaurea triumfetti subsp. aligera',
              'Hymenolobus procumbens var. revelieri',
              'Gnaphalium luteo-album',
              'Bromus hordeaceus subsp. thominii',
              'Carduncellus monspelliensium',
              'Centaurea triumfetti subsp. semidecurrens',
              'Odontites viscosa',
              'Ranunculus revelieri',
              'Odontites lutea',
              'Taraxacum braunblanquetii',
              'Galium pseudaristatum',
              'Odontites lanceolata',
              'Juncus acutus subsp. leopoldii')
}
sps_$matchType[sps_$initName%in%validated]="EXACT"
#sps_[sps_$matchType!="EXACT",]

sps_ = sps_[sps_$matchType=="EXACT",c('initName','species','speciesKey','confidence')]
colnames(sps_)=c('name_eva','species','speciesKey_GBIF','match_confidence')

setwd(eva)
write.table(sps_,'eva_species_list_match_gbif.csv',sep=";",row.names=F,col.names=T)

spCline=sp%>%
  merge(sps_[,c('name_eva','species')],by.x="Turboveg2.concept",by.y="name_eva",all.x=T)%>%
  select(PlotObservationID,species)%>%
  filter(!is.na(species))%>%
  unique()
PAeva=sites%>%
  filter(Location.uncertainty..m.<=100)%>%
  select(PlotObservationID,
         PlotID,
         TV2.relevé.number,
         Date.of.recording,
         Relevé.area..m..,
         Longitude,
         Latitude,
         Location.uncertainty..m.,
         Dataset,
         year,
         Access.regime)%>%
  merge(spCline,by="PlotObservationID",all.x=T)
colnames(PAeva)[c(3,4,5,8,9)]=c('TV2.releve.num',
                                'date',
                              'area.m2',
                                'coordinateUncertaintyInMeters',
                                'dataset')
setwd(main)
write.table(PAeva,'occ_PA_eva_CBNA_CBNmed.csv',sep=";",row.names=F,col.names=T)

# Extr CBN_a pour Sylvain
extCBNA=unique(PAeva[PAeva$year>=2016 & PAeva$dataset=="CBNA",colnames(PAeva)!="species"])
setwd(eva)
write.table(extCBNA,'extraction_eva_CBNA.csv',sep=";",row.names=F,col.names=T)

# Extr CBNmed pour Olivier
extCBNMed=unique(PAeva[PAeva$year>=2016 & PAeva$dataset=="CBNMed",colnames(PAeva)!="species"])
setwd(eva)
write.table(extCBNMed,'extraction_eva_CBNMed.csv',sep=";",row.names=F,col.names=T)


######
# Make PA_dispo table
######

setwd(eva)
cbn=read.csv('extraction_eva_CBNMed.csv',sep=";",header=T)
setwd(main)
PA=read.csv('occ3_sup2016_PA.csv',sep=";",header=T)

PA_ano=read.csv('PA_public_anonymised.csv',sep=";",header=T)

tmp=sample(1:dim(PA)[1],1)
cd=tmp:(tmp+40)

plot(PA$lon[cd],PA$lat[cd])
points(PA_ano$lon[cd],PA_ano$lat[cd],pch=4,col="red")


TMP=PA_ano%>%
  group_by(lon,lat,year,dayOfYear)%>%
  summarise()

######
# Data paper: Map of PO
######

require(ggmap)
po=read.csv(paste0(main,'Presences_only_train.csv'),sep=";",header=T)


mapo = get_stamenmap( bbox = c(left=-11,
                               bottom=32,
                               right=36,
                               top = 74), 
                      zoom = 5, maptype = "toner-lite")


scalo=4

toPlot=po%>%
  mutate(x=round(lon*scalo)/scalo,y=round(lat*scalo)/scalo)%>%
  group_by(x,y)%>%
  summarise(n=n())

breakos = c(-.1,quantile(as.vector(toPlot$n),prob=c(.05,.2,.4,.7,.8,.95)),max(toPlot$n))
breakos[2:length(breakos)]=round(breakos[2:length(breakos)])
toPlot$interval=cut(toPlot$n,breaks=breakos)
coloors = colorRampPalette(colors = c('goldenrod','darkorchid4'))(length(levels(toPlot$interval)))
leg = data.frame(num=levels(toPlot$interval),col=coloors)

p=ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=interval),alpha=.6)+
  scale_fill_manual(values=leg$col)+
  xlab('Longitude')+ylab('Latitude')+
  theme_bw()+
  theme(text=element_text(size=30))

setwd(main)
png(paste0('Figures/map_test.png'),height = 1000,width=1200)
print(p)
dev.off()

######
# Data paper: Map of PA
######


setwd(main)
test=read.csv('Presences_Absences_test.csv',sep=";",header=T)
blind=read.csv('test_blind.csv',sep=";",header=T)
train=read.csv('Presences_Absences_train.csv',sep=";",header=T)




gt=read.csv('gt.csv',sep=",",header=T)

all.equal(gt$Id,blind$Id)

cd=1:500
plot(blind$lon[cd],blind$lat[cd])

