# Shiny specific functions

# http://stackoverflow.com/questions/34530142/drop-down-checkbox-input-in-shiny
dropdownButton <- function(label = "", status = c("default"), ..., width = NULL) {
    status <- match.arg(status)
    # dropdown button content
    html_ul <- list(
        class = "dropdown-menu",
        style = if (!is.null(width))
            paste0("width: ", validateCssUnit(width), ";"),
        lapply(
            X = list(...),
            FUN = tags$li,
            style = "margin-left: 10px; margin-right: 10px;"
        )
    )
    # dropdown button apparence
    html_button <- list(
        class = paste0("btn btn-", status, " dropdown-toggle"),
        type = "button",
        `data-toggle` = "dropdown"
    )
    html_button <- c(html_button, list(label))
    html_button <- c(html_button, list(tags$span(class = "caret")))
    # final result
    tags$div(
        class = "dropdown",
        do.call(tags$button, html_button),
        do.call(tags$ul, html_ul),
        tags$script(
            "$('.dropdown-menu').click(function(e) {
            e.stopPropagation();
         });"))
}

#Custom input style
textInput3 <- function (inputId, label, value = "", ...){
    div(
    tags$label(label, `for` = inputId), 
    tags$input(id = inputId, type = "text", value = value, ...))
}


## List Parameters

uploadchoices <- list(
    "Loops/.rds" = 1,
    "ReadDepth/.bigWig" = 2,
    "ReadDepth/.bedgraph" = 3,
    "Methylation/.bigWig" = 4,
    "Methylation/.bedgraph" = 5,
    "Hi-C/.tgz" = 6
)

color.choices <- list(
    "Pallet" = 1,
    "Viridian" = 2,
    "Pewter" = 3,
    "Cerulean" = 4,
    "Vermillion" = 5,
    "Lavendar" = 6,
    "Celadon" = 7,
    "Fuchsia" = 8,
    "Saffron" = 9,
    "Cinnabar" = 10,
    "Indigo" = 11,
    "Master" = 12,
    "Black and Blue" = 13,
    "Heat" = 14,
    "Topology" = 15,
    "Blue to Red" = 16)


missingco.choices <- list(
    "Minimum" = "min",
    "Gray" = "gray97",
    "Black" = "black",
    "Red" = "red",
    "White" = "white")



## Useful Functions ##

.subsetRegion.quick <- function(loops, region, nanchors = 2) {
    g <- findOverlaps(region, loops@anchors)@to
    if(nanchors == 2) { cc <- loops@interactions[,1] %in% g & loops@interactions[,2] %in% g
    } else { cc <- xor(loops@interactions[,1] %in% g, loops@interactions[,2] %in% g) }
    
    ni <- matrix(loops@interactions[cc,], ncol = 2)
    colnames(ni) <- c("left", "right")
    nc <- matrix(loops@counts[cc,], ncol = dim(loops@counts)[2])
    colnames(nc) <- colnames(loops@counts)

    slot(loops, "interactions", check = TRUE) <- ni
    slot(loops, "counts", check = TRUE) <- nc
    slot(loops, "rowData", check = TRUE) <- loops@rowData[cc,, drop = FALSE]
    return(loops)
}
