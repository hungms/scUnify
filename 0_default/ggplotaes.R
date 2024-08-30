theme_border <- function(){
    theme(
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
        }

theme_line <- function(){
    list(
        scale_y_continuous(expand = c(0, 0)),
        theme(
            axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", size = 1)))
        }

facet_aes <- function(){
    theme(
	strip.background = element_blank(),
        strip.text = element_text(face="bold", size=12))}

umap_aes <- function(){
    list(
	xlab("UMAP1"),
        ylab("UMAP2"),
        theme(
            panel.background = element_rect(fill = "white"),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            aspect.ratio = 1))}
