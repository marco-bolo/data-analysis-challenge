#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
library(shinythemes)
library(bslib)
library(plyr)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(MetBrewer)
library(plotly)
library(patchwork)
library(pheatmap)
library(ggplotify)

# read in the data ####
# read in the metadata
metadata <- read.csv("app_data/submission_metadata_base_aquarium_selection.csv", sep=";",header=T)
row.names(metadata)=metadata$submissionID
head(metadata)

# the data frame for the PCA#submissionID the data frame for the PCA ####
pca_df <-read.csv("app_data/fish_pcoa_markerwise.csv",header=T)
head(pca_df)
# join the metadata and the pca
pca_plot<-join(pca_df,metadata,type="left",by=c("submissionID","submitterID","marker","referenceDB","pipeline"))


# the different data frames for the boxplots ####
bp_species <- read.csv("app_data/results_fish_targetClass_speciesLevel_alpha.csv", header=T)
bp_genus <- read.csv("app_data/results_fish_targetClass_genusLevel_alpha.csv", header=T)
bp_class_unit <- read.csv("app_data/results_fish_targetClass_unitLevel_alpha.csv", header=T)
bp_none_unit <- read.csv("app_data/results_fish_targetNone_unitLevel_alpha.csv", header=T)
# need to create a baseplot for the initial figure
bp_plot<-join(bp_species,metadata,type="left",by=c("submissionID","submitterID","marker","referenceDB","pipeline"), match="first")

# data input and processing for the heatmap ####
# read in the species guide
species_guide<-read.csv("app_data/aquarium_COI_12S_16S_allSpecies_genus.csv",header=T)
head(species_guide)
species.plot<-species_guide # create a new species guide for plotting
row.names(species.plot)=species.plot$assignment

# read in the data
raw.hm <- read_rds("app_data/megadata_no18s_250925.rds")

# make a list of the pipelines that found the target species
df.species <-filter(raw.hm, taxonomicLevel=="species" & assignment%in%species_guide$assignment)%>%
  group_by(submissionID,assignment)%>%
  mutate(uniqueID=paste(assignment,submissionID,sep="_"))
species_list<-unique(df.species$uniqueID)

# make a dataframe that has the data pipelines that found the target genus
df.genus <- filter(raw.hm, taxonomicLevel=="genus" & assignment%in%species_guide$genus) %>%
  group_by(submissionID,assignment)%>%
  mutate(uniqueID=paste(assignment,submissionID,sep="_"))
genus_list<-unique(df.genus$uniqueID)

# take the sum of the reps per pipeline and genus or species
df.filter<-filter(df.genus, (taxonomicLevel=="genus" & assignment%in%species_guide$genus) | (taxonomicLevel=="species" & assignment%in%species_guide$assignment)) %>%
  group_by(submissionID,taxonomicLevel,assignment)%>%
  summarise(readSum=sum(as.numeric(readAbundance), na.rm=T)) %>%
  mutate(unique_ID=paste(assignment,submissionID,sep="_"))

# now create a grid that includes all the species and the submissionID
# this creates the matrix including the possibility that a pipeline did not find a species
pipeline.df <-expand.grid(species_guide$assignment,unique(raw.hm$submissionID)) %>%
  as.data.frame()
names(pipeline.df)[1:2]=c("assignment","submissionID")
df.hm<-join(pipeline.df,select(species_guide,-Type),type="left") %>%
  mutate(uniqueID_specs=paste(assignment,submissionID,sep="_"), uniqueID_genus=paste(genus,submissionID,sep="_"),
         value=case_when(
           uniqueID_specs%in%species_list~2,
           uniqueID_genus%in%genus_list~1,
           TRUE~0
         )
  )

head(df.hm)
df.hm.plot <- dcast(df.hm, submissionID~assignment,value.var = "value")
#add the submissionIDs as rownames
rownames(df.hm.plot)=df.hm.plot$submissionID

# Define UI for application####
ui <- fluidPage(
  theme=shinytheme("cerulean"),
  page_navbar(
  # Application title
  title=div(img(src="MARCO-BOLO_logo_col.jpg",height = "57.5px", width = "auto"),"WP2.2- Data Analysis Challenge"),
  tabPanel("Fish-PCoA",
           sidebarLayout(
             sidebarPanel(
               width=2,
               fluid=TRUE,
               h4("PCOA"), # title of the next section
               # choose the variable that is color of the PCoA based on the df
               varSelectInput("columnVar", "Color by:", pca_plot, selected = "submitterID")
             ),
             mainPanel(
               plotlyOutput("plot",width = "80vw", height = "85vh")
             )
           )
  ),
  tabPanel("Fish-Barplot of diversity metrics",
           sidebarLayout(
             sidebarPanel(
               width=2,
               h4("Diversity Boxplot"),  # title of the section, for boxplot
               # choosing the dataset used
               selectInput("bp_df","Taxonomy Level:",c("Target Species Level","Target Genus Level","Target Unit Level","Non-target Unit Level"), selected="Target Species Level"),
               # choosing the x variable for boxplot
               varSelectInput("xvar", "X variable", bp_plot, selected = "submitterID"),
               # choosing the y variable for boxplot
               varSelectInput("yvar", "Y variable", bp_plot, selected = "shannon"),
               # choosing the color for boxplot  by the df
               varSelectInput("columnVar2", "Color by:", bp_plot, selected = "submitterID")
             ),
             mainPanel(
               plotlyOutput("plot2",width = "80vw", height = "85vh")
             )
           )
  ),
  tabPanel("Fish-Pipeline accuracy heat map",
           sidebarLayout(
             sidebarPanel(
               width=2,
               h4("Heatmap of species found by pipeline"),     # title
               # choosing the dataset used
               selectInput("hm_df","Taxonomy Level:",c("Species Level","Genus Level"), selected="Species Level"),
               # select grouping metric
               selectInput("columnVar3","Clustering metric:",c("submissionID","submitterID","pipeline","referenceDB","referenceDB_simplified",
                                                               "ASV.OTU","screen_other","long_formatter","tax_assign_cat"), selected="ASV.OTU")
             ),
             mainPanel(
               plotOutput("plot3",width = "80vw", height = "160vh")
             )
           )
  )
)
)
# Define server logic  ####
server <- function(input, output) {
  
    output$plot <- renderPlotly({
    # pcoa figures
    theme_ordination <- function() {
      list(theme_classic() + # base
             theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                   axis.title = element_text(size = 12, face = "bold"),
                   axis.text = element_text(size = 10),
                   axis.line = element_line(color = "black"),
                   panel.grid = element_blank()),
           # dashed zero lines
           geom_hline(yintercept = 0, linetype = "dashed", color = "grey50"),
           geom_vline(xintercept = 0, linetype = "dashed", color = "grey50"))}
    
    p <- ggplot(pca_plot, aes(x = PC1, y = PC2, color = !!input$columnVar, label=submissionID)) +
      geom_point(size = 3) +
      facet_wrap(~marker) +
      theme_ordination() +
      labs(title = "PCoA (Bray-Curtis)", x = "PC1", y = "PC2") +
      scale_color_discrete(guide = guide_legend(ncol = 1)) # +
    
    ggplotly(p) 
  })
    # select the dataframe for the boxplot
    selected_data = reactive({
      switch(input$bp_df,
             "Target Species Level"= bp_species,
             "Target Genus Level"= bp_genus,
             "Target Unit Level" = bp_class_unit,
             "Non-target Unit Level" = bp_none_unit
      )
    })
    
    # make a plotly boxplot
    output$plot2 <- renderPlotly({
      # use reactive expression selected_data() to get the chosen data frame
      bp_plot<-join(selected_data(),metadata,type="left",by=c("submissionID","submitterID","marker","referenceDB","pipeline"))
      
      # create ggplot object
      p<-ggplot(data = bp_plot, aes(y= !!input$yvar, x = !!input$xvar, color = !!input$columnVar2)) +
        geom_boxplot(position = "dodge2") +
        geom_point(position = position_jitterdodge(),size=.25, alpha=0.8) +
        facet_grid(marker ~.) +
        ggtitle("Species richness (for assignments to target class)") +
        scale_color_discrete(guide = guide_legend(ncol = 1)) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.title = element_text(face = "bold"),
              #axis.text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.line = element_line(color = "black"),
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = 30))
      
      ggplotly(p) 
      
    })
  
    output$plot3 <- renderPlot({
      if (input$hm_df=="Species Level"){
        df.hm.plot <- dcast(df.hm, submissionID~assignment,value.var = "value") # recast the dataframe by species
        species.plot<-species_guide # create a new species guide for plotting
        row.names(species.plot)=species.plot$assignment
        
      }
      if (input$hm_df=="Genus Level"){
        df.hm.plot <- dcast(df.hm, submissionID~genus,value.var = "value",fun.aggregate=max) # recast the dataframe by genus
        species.plot<-unique(species_guide[c(1,3)]) # create a new sepcies guide for plotting that reduces the redundancy of the genus
        row.names(species.plot)=species.plot$genus
        
      }
      #add the submissionIDs as rownames
      rownames(df.hm.plot)=df.hm.plot$submissionID
    
      # heatmap plotting function, input is the data frame, metadata
      # what variable is to be clustered by
      hmPlot <- function(df1,metadata,species.plot,cluster){
        markers_list=c("12S","16S","COI") # make a list of the markers
        plot_markers = list() # create an empty list for plotting
        i = 1
        for (mark in markers_list){
          sub <- filter(df1, (submissionID%in%filter(metadata, marker==mark)$submissionID))
          sub.plot = sub[-1]
          p <- pheatmap(sub.plot,color=colorRampPalette(brewer.pal(n = 3, name ="Greens"))(3), 
                        annotation_row = metadata[cluster],annotation_col=species.plot[1],
                        fontsize_col=8,cellheight=8,fontsize_row=8,
                        legend_breaks=c(0,1,2),legend_labels=c("Not Detected","Genus Detected","Species Detected"),
                        angle_col="45",
                        main=mark)
          
          plot_markers[[i]]=p[[4]]
          i = i+1
        }
        return(plot_markers)
      }
     plot_list<-hmPlot(df.hm.plot,metadata,species.plot,input$columnVar3)
     as.ggplot(plot_list[[1]]) / as.ggplot(plot_list[[2]]) / as.ggplot(plot_list[[3]]) + plot_layout(width=c())
     })
}

# Run the application ####
shinyApp(ui = ui, server = server)