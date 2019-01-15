library(rgeos)
library(sp)
library(raster)

#####
# Functions
#####

# Fonction: Taxon name -> TaxrefName
getExactMatchTaxref = function(scName,ref,rangs){
  scName = toupper(scName)
  
  # matching dans ref
  id = which(ref$NOM_COMPLET==scName)
  # Y a t'il eu un matching sur Nom+auteur+année?
  if(length(id)==0){
    #print('no match with author')
    # Non, on tente avec Nom + auteur
    id = which(ref$NOM_AUTEUR==scName)
    if(length(id)==0){
      # Non, on tente Nom seul
      nameAlone = paste(strsplit(as.character(scName)," ")[[1]][1:2],collapse = " ")
      id = which(ref$LB_NOM==nameAlone)
    }
  }
  # n'y a t'il eu aucun nom matché?
  if(length(id)==0){
    # non aucun: le nom est inconnu du référentiel
    #print(paste('No match at all for',scName))
    valid_name = NA
    # sinon, on détermine combien de 
    # noms valides d'espèces correspondent
    # aux noms matchés
  }else{
    valid_match = NULL
    for(ID in id){
      niveau = rangs$RG_LEVEL[rangs$RANG==ref$RANG[ID]]
      if(niveau>=290){
        id_es = reach.higher.rank.TAXREF(ID,ref,rangs)
        if(!is.na(id_es)){
          valid_match = c(valid_match,ref$NOM_VALIDE[id_es])
        }
      }
      # Si tous les match de niveau supérieur
      # ou égal à l'espèce mènent à la MEME espèce de référence
      # on prend celle là
      # Sinon, c'est qu'il y a une ambiguité, on laisse tomber
      if(length(unique(valid_match))==1){
        valid_name = unique(valid_match)
      }else{valid_name=NA}
    }
  }
  
  return(valid_name)}

# Function: index taxon -> index taxon corresponding species
reach.higher.rank.TAXREF = function(id,ref,rangs,code=290){
  niveau = rangs$RG_LEVEL[rangs$RANG==ref$RANG[id]]
  if(niveau==code){
    return(id)
  }else if(niveau<code){
    return(NA)
  }else{
    # On va chercher le taxon de référence correspondant à celui là
    id_ref = which(ref$CD_NOM == ref$CD_REF[id])
    # On va chercher le taxon parent du taxon de ce taxon de référence
    id_sup = which(ref$CD_NOM == ref$CD_SUP[id_ref]) 
    return(reach.higher.rank.TAXREF(id_sup,ref,rangs,code))
  }
}

# Function: occurrences -> filtered occurrences over French metropolitan territory
filter_FR = function(occ,M,proj_ini = "+proj=longlat +datum=WGS84",project=F){
  #Load bordering countries (use online connection)
  border_countries = c( "FRA", "ESP", "ITA", "BEL" , "LUX", "CHE", "DEU" )
  for (pays in border_countries){
    if (pays=="FRA"){ MAPS=   spTransform( getData('GADM', country=pays, level=0),proj ) 
    }else{ MAPS = rbind(MAPS, spTransform( getData('GADM', country=pays, level=0),proj ) ) }
  }
  france=MAPS[1,]
  
  # We project WGS84 occurrences to LCC for measuring metric distances between
  # occurrences location and France polygon 
  coos = data.frame(Longitude=occ$Longitude,Latitude=occ$Latitude)
  pts_occ = SpatialPointsDataFrame( coords=coos,data=data.frame(id = 1:dim(occ)[1]),proj4string = CRS(proj_ini))
  pts_occ = spTransform(pts_occ,CRS("+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"))

  # distance to France polygon
  nBand = 1000
  n = dim(pts_occ)[1]
  band = n %/% nBand;rest = n %% nBand;
  indexes = lapply(1:nBand , function(i) (i-1)*band + (1:band) )
  indexes[[length(indexes)+1]] = length(indexes)*band + (1:rest)
  dist_fr = NULL
  for(i in 1:length(indexes)){
    dist_fr = c(dist_fr , gDistance(france, pts_occ[indexes[[i]],],byid=T)[,1] )
    saveRDS(dist_fr,'C:/Users/Christophe/pCloud local/0_These/data/geolifeclef/GLC2019/occurrences/dist_fr')
    flush.console()
    cat('   \r    ',100*i/length(indexes),'%     \r')
  }
  
  # identify points at distance between 0 and M meters from the border
  not_in_but_close = dist_fr<=M & dist_fr>0
  
  # For those, identify those closer from France than other neighboor country
  autres_pays=MAPS[2:length(MAPS),]
  if(sum(not_in_but_close)>0 & project){
    dist_others = gDistance(autres_pays,pts_occ[not_in_but_close,],byid=T)
    dist_others = sapply( 1:dim(dist_others)[1],function(i) min(dist_others[i,]) )
    to_project = rep(F,length(dist_fr))
    to_project[not_in_but_close] = dist_others > dist_fr[not_in_but_close]
    tokeep = dist_fr==0 | to_project
  }else{
    tokeep = dist_fr<=M
  }
  
  new_occ= occ
  new_occ = new_occ[tokeep,]
  return(new_occ)
}


#####
# Preparing final train Datasets
# (i) Taxonomic matching with Taxref
# (ii) spatial filtering FR
# (iii) shuffling 
#####

# Directory where to read unmatched occurrences datasets
dir = "C:/Users/Christophe/pCloud local/0_These/data/geolifeclef/GLC2019/occurrences/"
# Each Source dataset is a ";" separated .csv file with headers containing 
# at least 3 columns :
# "Longitude" : Longitude coordinate WGS84
# "Latitude": Latitude coordinate WGS84 
# "scName" : scientific name of the observed taxon  
Sources = c('gbif_2018_0.csv',
            'PL_complete_0.csv',
            'GLC_2018_0.csv')

taxaName_glc19SpId = data.frame(taxaName=NA,glc19SpId=NA)
taxaName_glc19SpId =   taxaName_glc19SpId[-1,,drop=F]

# load Taxref 12 referential
TaxrefDir = "C:/Users/Christophe/pCloud local/0_These/data/Taxonomy/référentiels/Taxref12.0/"
ref = read.csv(paste(TaxrefDir,'TAXREFv12_Plantae.csv',sep=""),sep=";",header=T)
ref$NOM_COMPLET = as.character(ref$NOM_COMPLET)
ref$NOM_COMPLET = toupper(ref$NOM_COMPLET)
ref$LB_NOM = as.character(ref$LB_NOM)
ref$LB_NOM = toupper(ref$LB_NOM)
ref$NOM_VALIDE = as.character(ref$NOM_VALIDE)
ref$RANG = as.character(ref$RANG)
# We isolate Author name from year
tmp = strsplit(as.character(ref$LB_AUTEUR),",")
authorsAlone = sapply(tmp,function(v) v[1])
ref$AUTEUR = toupper(authorsAlone)
ref$NOM_AUTEUR = toupper(paste(ref$LB_NOM,ref$AUTEUR))

# Tableau des RANG 
rangs = read.csv(paste(TaxrefDir,'rangs_note.csv',sep=""),sep=";",header=T)
rangs$RANG=as.character(rangs$RANG)
rangs$RG_LEVEL=as.numeric(rangs$RG_LEVEL)

# Prepare geographic filtering 
M=30
proj_ini = "+proj=longlat +datum=WGS84"
project=F

for(i in 1:length(Sources)){
  print(Sources[i])
  tab = read.csv(paste(dir,Sources[i],sep=""),sep=";",header=T)
  tab$glc19SpId = NA
  
  scNames = unique(tab$scName)
  # Make Names table (if not already on disk)
  setwd(dir)
  if(!'scName_TaxrefName.csv'%in%list.files()){
    scName_TaxrefName = data.frame(scName=NA,TaxrefName=NA)
    scName_TaxrefName = scName_TaxrefName[-1,,drop=F]
  }else{
    scName_TaxrefName = read.csv("scName_TaxrefName.csv",sep=";",header=T)
  }
  
  # Taxonomic matching with Taxref, and attribution of glc19SpId
  for(name in scNames){
    k = which(scNames==name)
    # On vérifie d'abord si le taxon
    # n'a pas déjà été traité dans scName_TaxrefName
    if(name%in%scName_TaxrefName$scName){
      cd = scName_TaxrefName$scName==name
      TaxrefName = scName_TaxrefName$TaxrefName[cd]
      if(!is.na(TaxrefName)){
        cd = taxaName_glc19SpId$taxaName==TaxrefName
        tab$glc19SpId[tab$scName==name] = taxaName_glc19SpId$glc19SpId[cd]
      }
    }else{
      TaxrefName = getExactMatchTaxref(name,ref,rangs)
      scName_TaxrefName = rbind( scName_TaxrefName, data.frame(scName=name,TaxrefName=TaxrefName) )
      if(!is.na(TaxrefName) & !TaxrefName%in%taxaName_glc19SpId$taxaName){
        newId = dim(taxaName_glc19SpId)[1]+1
        taxaName_glc19SpId = rbind(taxaName_glc19SpId , data.frame(taxaName=TaxrefName,glc19SpId=newId,test=F) )
      }else if(!is.na(TaxrefName)){
        newId = taxaName_glc19SpId$glc19SpId[taxaName_glc19SpId$taxaName==TaxrefName]
      }else{
        newId = NA
      }
      tab$glc19SpId[tab$scName==name]=newId
    }
    if(k/100 == round(k/100)){
      flush.console()
      cat('    \r  ',100*k/length(scNames),'%       \r   ')
    }
  }
  setwd(dir)
  write.table(scName_TaxrefName,'scName_TaxrefName.csv',sep=";",row.names=F,col.names=T)
  
  nTaxref = length(unique(tab$glc19SpId[!is.na(tab$glc19SpId)]))
  print(paste('Taxon recovery rate in Taxref',nTaxref/length(scNames)))
  
  tabo = tab[!is.na(tab$glc19SpId),,drop=F]
  
  # French borders filter
  taboo = filter_FR(tabo,M)
  
  # Dataset shuffling 
  taboo = taboo[sample(1:dim(taboo)[1],dim(taboo)[1]),,drop=F]
  
  # Save post-treated dataset
  write.table(taboo,
              paste(dir,substr(Sources[i],1,nchar(Sources[i])-6),'.csv',sep=""),
              sep=";",col.names=T,row.names=F)
}
setwd(dir)
write.table(taxaName_glc19SpId,'taxaName_glc19SpId.csv',sep=";",row.names=F,col.names=T)



