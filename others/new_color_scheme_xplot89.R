
xplot89 <- function(muts_basis,text_size,plot_title){
  
  indel_type_4_figurelabel <-
    structure(
      list(
        IndelType =
          c(
            
            
            
            ## === SINGLE BASE DEL ===
            # Single base deletions -  C
            "[Del(C):R1]A", "[Del(C):R1]T", "[Del(C):R2]A", "[Del(C):R2]T",
            "[Del(C):R3]A", "[Del(C):R3]T", "[Del(C):R(4,5)]A", "[Del(C):R(4,5)]T",
            "[Del(C):R(1,5)]G", "Del(C):R(6,9)",
            
            # Single base deletions - T
            "A[Del(T):R(1,4)]A", "A[Del(T):R(1,4)]C",
            "A[Del(T):R(1,4)]G", "C[Del(T):R(1,4)]A", "C[Del(T):R(1,4)]C",
            "C[Del(T):R(1,4)]G", "G[Del(T):R(1,4)]A", "G[Del(T):R(1,4)]C",
            "G[Del(T):R(1,4)]G", "A[Del(T):R(5,7)]A", "A[Del(T):R(5,7)]C",
            "A[Del(T):R(5,7)]G", "C[Del(T):R(5,7)]A", "C[Del(T):R(5,7)]C",
            "C[Del(T):R(5,7)]G", "G[Del(T):R(5,7)]A", "G[Del(T):R(5,7)]C",
            "G[Del(T):R(5,7)]G", "A[Del(T):R(8,9)]A", "A[Del(T):R(8,9)]C",
            "A[Del(T):R(8,9)]G", "C[Del(T):R(8,9)]A", "C[Del(T):R(8,9)]C",
            "C[Del(T):R(8,9)]G", "G[Del(T):R(8,9)]A", "G[Del(T):R(8,9)]C",
            "G[Del(T):R(8,9)]G",
            ## === END SINGLE BASE DEL ===
            
            
            
            # === Single base INS ===
            # Single base insertions - C (5)
            "A[Ins(C):R0]A",
            "A[Ins(C):R0]T",
            "Ins(C):R(0,3)",
            "Ins(C):R(4,6)",
            "Ins(C):R(7,9)",
            
            # Single base insertions - T (24)
            "A[Ins(T):R(0,4)]A",
            "A[Ins(T):R(0,4)]C", "A[Ins(T):R(0,4)]G", "C[Ins(T):R(0,4)]A",
            "C[Ins(T):R(0,4)]C", "C[Ins(T):R(0,4)]G", "G[Ins(T):R(0,4)]A",
            "G[Ins(T):R(0,4)]C", "G[Ins(T):R(0,4)]G", "A[Ins(T):R(5,7)]A",
            "A[Ins(T):R(5,7)]C", "A[Ins(T):R(5,7)]G", "C[Ins(T):R(5,7)]A",
            "C[Ins(T):R(5,7)]C", "C[Ins(T):R(5,7)]G", "G[Ins(T):R(5,7)]A",
            "G[Ins(T):R(5,7)]C", "G[Ins(T):R(5,7)]G", "A[Ins(T):R(8,9)]A",
            "A[Ins(T):R(8,9)]C", "A[Ins(T):R(8,9)]G", "C[Ins(T):R(8,9)]A",
            "C[Ins(T):R(8,9)]C", "C[Ins(T):R(8,9)]G", "G[Ins(T):R(8,9)]A",
            "G[Ins(T):R(8,9)]C", "G[Ins(T):R(8,9)]G",
            # === END Single base INS ===
            
            
            # Longer deletions no MH (6)
            "Del(2,4):R1",
            "Del(5,):R1",
            "Del(2,8):U(1,2):R(2,4)",
            "Del(2,):U(1,2):R(5,9)",
            "Del(3,):U(3,):R2",
            "Del(3,):U(3,):R(3,9)",
            
            
            
            # Longer insertions (6)
            "Ins(2,4):R0", "Ins(5,):R0",
            "Ins(2,4):R1", "Ins(5,):R1",
            "Ins(2,):R(2,4)", "Ins(2,):R(5,9)",
            
            
            
            # Deletions with MH (7)
            "Del(2,5):M1", "Del(3,5):M2", "Del(4,5):M(3,4)", "Del(6,):M1",
            "Del(6,):M2", "Del(6,):M3", "Del(6,):M(4,)",
            
            "Complex"),
        
        Indel =
          c(
            
            rep("Del(C)", 10),
            rep("Del(T)", 27),
            
            
            rep("Ins(C)", 5),
            rep("Ins(T)", 27),
            
            
            # Six long deletions
            rep("Del(2,):R(0,9)", 6),
            
            # Six long insertions
            rep("Ins(2,)", 6),
            
            
            
            # Seven deletions with MH
            rep("Del(2,):M(1,)", 7),
            
            "Complex"),
        
        
        Indel3 =
          c(
            
            
            # 1 bp deletions, 10 + 27
            rep("Deletion", 10 + 27),
            
            # 1 bp insertions, 5 + 27
            rep("Insertion", 5 + 27),
            
            
            
            # 6 long del, no MH
            rep("Deletion", 6),
            
            # Six long insertions
            rep("Insertion", 6),
            
            
            
            # 7 long del, MH
            rep("Deletion", 7),
            
            "Complex"),
        
        Figlabel =
          c(
            
            ## === 1 BP DELETION
            # Del C
            "[C1]A", "[C1]T", "[C2]A", "[C2]T", "[C3]A", "[C3]T", "[C(4,5)]A",
            "[C(4,5)]T", "[C(1,5)]G", "C(6,9)",
            
            # Del T
            "A[T(1,4)]A", "A[T(1,4)]C",
            "A[T(1,4)]G", "C[T(1,4)]A", "C[T(1,4)]C", "C[T(1,4)]G", "G[T(1,4)]A",
            "G[T(1,4)]C", "G[T(1,4)]G", "A[T(5,7)]A", "A[T(5,7)]C", "A[T(5,7)]G",
            "C[T(5,7)]A", "C[T(5,7)]C", "C[T(5,7)]G", "G[T(5,7)]A", "G[T(5,7)]C",
            "G[T(5,7)]G", "A[T(8,9)]A", "A[T(8,9)]C", "A[T(8,9)]G", "C[T(8,9)]A",
            "C[T(8,9)]C", "C[T(8,9)]G", "G[T(8,9)]A", "G[T(8,9)]C", "G[T(8,9)]G",
            
            ## === END 1 BP DELETION
            
            
            # Ins C
            
            "A[C0]A", "A[C0]T",
            "C(0,3)", "C(4,6)", "C(7,9)",
            
            # Ins T
            
            "A[T(0,4)]A", "A[T(0,4)]C", "A[T(0,4)]G",
            "C[T(0,4)]A", "C[T(0,4)]C", "C[T(0,4)]G", "G[T(0,4)]A", "G[T(0,4)]C",
            "G[T(0,4)]G", "A[T(5,7)]A", "A[T(5,7)]C", "A[T(5,7)[G", "C[T(5,7)[A",
            "C[T(5,7)[C", "C[T(5,7)[G", "G[T(5,7)[A", "G[T(5,7)[C", "G[T(5,7)]G",
            "A[T(8,9)]A", "A[T(8,9)]C", "A[T(8,9)]G", "C[T(8,9)]A", "C[T(8,9)]C",
            "C[T(8,9)]G", "G[T(8,9)]A", "G[T(8,9)]C", "G[T(8,9)]G",
            
            # Long del no MH
            "L(2,4):R1", "L(5, ):R1",
            "L(2,8):U(1,2):R(2,4)", "L(2, ):U(1,2):R(5,9)",
            "L(3, ):U(3,):R2", "L(3, ):U(3,):R(3,9)",
            
            
            # Long insertion
            "L(2,4):R0",
            "L(5, ):R0",
            "L(2,4):R1",
            "L(5, ):R1",
            "L(2, ):R(2,4)",
            "L(2, ):R(5,9)",
            
            
            
            # Long del MH
            "L(2,5):M1", "L(3,5):M2",
            "L(4,5):M(3,4)", "L(6, ):M1",
            "L(6, ):M2", "L(6, ):M3",
            "L(6, ):M(4, )",
            
            
            "Complex")),
      
      
      class = "data.frame", row.names = c(NA, -89L))
  
  muts_basis_melt <- reshape2::melt(muts_basis,"IndelType")
  
  muts_basis_melt <- merge(indel_type_4_figurelabel, muts_basis_melt,by="IndelType",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("IndelType","Indel","Indel3","Figlabel","Sample","freq")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)
  
  indel_mypalette_fill <- c("#000000", # FEABB9 Complex
                            "#61409b", # Del(2,):M(1,)
                            "#f14432", # Del(2,):R(0,9) #EE6677
                            "#fdbe6f", # Del(C)
                            "#ff8001", # Del(T)
                            "#4a98c9", #"#668D3C"  Ins(2,)
                            "#b0dd8b", #"#007996"  Ins(C)
                            "#36a12e") # Ins(T)
  
  
  
  
  indel_positions <- indel_type_4_figurelabel$IndelType
  indel_positions_labels <- indel_type_4_figurelabel$Figlabel
  
  # color blocks for indel bases
  
  entry <- table(indel_type_4_figurelabel$Indel)
  
  
  order_entry <- c(
    "Del(C)", "Del(T)",
    "Ins(C)", "Ins(T)",
    "Del(2,):R(0,9)",
    "Ins(2,)",
    "Del(2,):M(1,)",
    "Complex")
  entry <- entry[order_entry]
  blocks <- data.frame(Type=unique(indel_type_4_figurelabel$Indel),
                       fill=indel_mypalette_fill,
                       xmin=c(0,cumsum(entry)[-length(entry)])+0.5,
                       xmax=cumsum(entry)+0.5)
  blocks$ymin <- max(muts_basis_melt$freq)*1.08#
  blocks$ymax <- max(muts_basis_melt$freq)*1.2
  # blocks$labels <-c("1bp C", "1bp T", ">=2bp", "1bp C", "1bp T", ">=2bp", "Mh", "X")
  blocks$labels = c("1bp C", "1bp T", "1bp C", "1bp T", ">=2bp", ">=2bp", "Mh", "X")
  blocks$cl <-c("black", "black", "black", "black", "white", "white", "white",  "white")
  
  # grey blocks for insertion or deletion or complex
  indel_mypalette_fill3 <- c("#000000", # FEABB9 Complex
                             "#888888", # Deletion
                             "#DDDDDD") # Insertion
  
  
  blocks3 = blocks
  blocks3 %>%
    dplyr::mutate(Type = ifelse(Type %in% c("Del(C)", "Del(T)"), "Del1", Type)) %>%
    dplyr::mutate(Type = ifelse(Type %in% c("Ins(C)", "Ins(T)"), "Ins1", Type)) %>%
    dplyr::group_by(Type) %>%
    dplyr::summarise(xmin = min(xmin), xmax = max(xmax)) %>%
    dplyr::mutate(labels = substr(Type, 1, 3)) %>%
    dplyr::mutate(labels = ifelse(labels == "Com", "X", labels)) %>%
    dplyr::arrange(xmin) -> blocks3
  blocks3$ymin <- max(muts_basis_melt$freq)*1.2
  blocks3$ymax <- max(muts_basis_melt$freq)*1.32
  blocks3$cl <- "white"
  #blocks3$Type = blocks3$labels
  blocks3$Type =c("Del1","Ins1","Del2","Ins2","DelMH","X")
  blocks3$cl[1:2] <- "black"
  indel_mypalette_fill_all <- c("Del1" = "#fe9f38",
                                "Ins1" = "#73bf5d",
                                "Del2" = "#f14432",
                                "Ins2" = "#4a98c9",
                                "DelMH" = "#61409b",
                                "X" = "black",
                                "Del(2,):M(1,)"= "#61409b", # Deletion
                                "Del(2,):R(0,9)" = "#f14432", # Del(2,):R(0,9) #EE6677
                                "Del(C)" = "#fdbe6f", # Del(C)
                                "Del(T)" = "#ff8001", # Del(T)
                                
                                "Ins(2,)" = "#4a98c9", #"#668D3C"  Ins(2,)
                                "Ins(C)"="#b0dd8b", #"#007996"  Ins(C)
                                "Ins(T)" = "#36a12e",
                                "Complex" = "black"# Ins(T)
  ) # Insertion
  
  
  p <- ggplot2::ggplot(data=muts_basis_melt, ggplot2::aes(x=IndelType, y=freq,fill=Indel)) +
    ggplot2::geom_bar(stat="identity",position="dodge", width=.7)+ggplot2::xlab("Indel Types") +
    ggplot2::ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  
  p <- p+ggplot2::scale_x_discrete(limits = indel_positions,labels = indel_positions_labels) +
    ggplot2::ggtitle(plot_title)
  
  p <- p+ggplot2::scale_fill_manual(values= indel_mypalette_fill_all)+
    ggplot2::scale_y_continuous(limits = c(0,unique(blocks3$ymax)),
                                labels=scales::number_format(accuracy = 0.01), expand = c(0,0))
  
  p <- p+ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, size=5,colour = "black",hjust=1),
      axis.text.y=ggplot2::element_text(size=10,colour = "black"),
      #    axis.line.y=element_blank(),
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size=15),
      axis.title.y = ggplot2::element_text(size=15))
  
  ## Add the overhead blocks for insertion or deletion or complex
  p <- p+ggplot2::geom_rect(
    data = blocks3,
    ggplot2::aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=Type,colour = "white"),
    inherit.aes = F) +
    # geom_text(data=blocks,aes(x=(xmax+xmin)/2,y=(ymax+ymin)/2,label=labels),size=text_size,inherit.aes = F,colour="white")
    ggplot2::geom_text(data=blocks3,
                       ggplot2::aes(x=(xmax+xmin)/2,y=(ymax+ymin)/2,label=labels, colour=cl),
                       size=text_size,fontface="bold",inherit.aes = F)+
    ggplot2::scale_colour_manual(values=c("black", "white"))
  
  ## Add the overhead blocks for indel bases
  p <- p+ggplot2::geom_rect(data = blocks,
                            ggplot2::aes(xmin=xmin,ymin=ymin,xmax=xmax,
                                         ymax=ymax,fill=Type,colour = "white"),inherit.aes = F)+
    ggplot2::geom_text(data=blocks,
                       ggplot2::aes(x=(xmax+xmin)/2,y=(ymax+ymin)/2,label=labels, colour=cl),
                       size=text_size,fontface="bold",inherit.aes = F)
  
  return(p)
  
}