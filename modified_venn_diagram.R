venn.diagram <- function (x, filename, height = 3000, width = 3000, resolution = 500, 
                          imagetype = "tiff", units = "px", compression = "lzw", na = "stop", 
                          main = NULL, sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
                          main.fontfamily = "serif", main.col = "black", main.cex = 1, 
                          main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                          sub.fontfamily = "serif", sub.col = "black", sub.cex = 0.85, 
                          sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE, 
                          print.mode = "raw", sigdigs = 3, big.mark = ",", direct.area = FALSE, area.vector = 0, 
                          hyper.test = FALSE, total.population = NULL, lower.tail = TRUE, 
                          ...) 
{
  time.string = gsub(":", "-", gsub(" ", "_", as.character(Sys.time())))
  if (!is.null(filename)) {
    flog.appender(appender.file(paste0(filename, ".", time.string, 
                                       ".log")), name = "VennDiagramLogger")
  }
  else {
    flog.appender(appender.file(paste0("VennDiagram", time.string, 
                                       ".log")), name = "VennDiagramLogger")
  }
  out.list = as.list(sys.call())
  out.list[[1]] <- NULL
  out.string = capture.output(out.list)
  flog.info(out.string, name = "VennDiagramLogger")
  if (direct.area) {
    if (1 == length(area.vector)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = area.vector[1], 
                                                 category = list.names, ind = FALSE, ...)
    }
    if (3 == length(area.vector)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = area.vector[1], 
                                                   area2 = area.vector[2], cross.area = area.vector[3], 
                                                   category = category.names, ind = FALSE, print.mode = print.mode, 
                                                   sigdigs = sigdigs, big.mark = ",", ...)
    }
    if (7 == length(area.vector)) {
      grob.list <- VennDiagram::draw.triple.venn(area1 = 0, 
                                                 area2 = 0, area3 = 0, n12 = 0, n23 = 0, n13 = 0, 
                                                 n123 = 0, category = category.names, ind = FALSE, 
                                                 list.order = 1:3, print.mode = print.mode, sigdigs = sigdigs, big.mark = ",",
                                                 area.vector = area.vector, direct.area = TRUE, 
                                                 ...)
    }
    if (15 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quad.venn(area1 = 0, 
                                               area2 = 0, area3 = 0, area4 = 0, n12 = 0, n13 = 0, 
                                               n14 = 0, n23 = 0, n24 = 0, n34 = 0, n123 = 0, 
                                               n124 = 0, n134 = 0, n234 = 0, n1234 = 0, category = category.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, big.mark = ",",
                                               area.vector = area.vector, direct.area = TRUE, 
                                               ...)
    }
    if (31 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = 0, 
                                                    area2 = 0, area3 = 0, area4 = 0, area5 = 0, n12 = 0, 
                                                    n13 = 0, n14 = 0, n15 = 0, n23 = 0, n24 = 0, 
                                                    n25 = 0, n34 = 0, n35 = 0, n45 = 0, n123 = 0, 
                                                    n124 = 0, n125 = 0, n134 = 0, n135 = 0, n145 = 0, 
                                                    n234 = 0, n235 = 0, n245 = 0, n345 = 0, n1234 = 0, 
                                                    n1235 = 0, n1245 = 0, n1345 = 0, n2345 = 0, n12345 = 0, 
                                                    category = category.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, big.mark = ",", area.vector = area.vector, 
                                                    direct.area = TRUE, ...)
      format(grob.list, big.mark = ",", scientific = FALSE)
    }
  }
  else {
    if (force.unique) {
      for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
      }
    }
    if ("none" == na) {
      x <- x
    }
    else if ("stop" == na) {
      for (i in 1:length(x)) {
        if (any(is.na(x[[i]]))) {
          flog.error("NAs in dataset", call. = FALSE, 
                     name = "VennDiagramLogger")
          stop("NAs in dataset", call. = FALSE)
        }
      }
    }
    else if ("remove" == na) {
      for (i in 1:length(x)) {
        x[[i]] <- x[[i]][!is.na(x[[i]])]
      }
    }
    else {
      flog.error("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"", 
                 name = "VennDiagramLogger")
      stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
    }
    if (0 == length(x) | length(x) > 5) {
      flog.error("Incorrect number of elements.", call. = FALSE, 
                 name = "VennDiagramLogger")
      stop("Incorrect number of elements.", call. = FALSE)
    }
    if (1 == length(x)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
                                                 category = list.names, ind = FALSE, ...)
    }
    else if (2 == length(x)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
                                                   area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
                                                                                                         x[[2]])), category = category.names, ind = FALSE, 
                                                   print.mode = print.mode, sigdigs = sigdigs, big.mark = ",", ...)
      idx <- sapply(p, function(i) grepl("text", i$name))
      
      for(i in 1:3){
        p[idx][[i]]$label <- 
          format(as.numeric(p[idx][[i]]$label), big.mark=",", scientific=FALSE)
      }
      
      grid.newpage()
      grid.draw(p) 
    }
    else if (3 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      list.names <- category.names
      nab <- intersect(A, B)
      nbc <- intersect(B, C)
      nac <- intersect(A, C)
      nabc <- intersect(nab, C)
      grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
                                                 area2 = length(B), area3 = length(C), n12 = length(nab), 
                                                 n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
                                                 category = list.names, ind = FALSE, list.order = 1:3, 
                                                 print.mode = print.mode, sigdigs = sigdigs, big.mark = ",", ...)
    }
    else if (4 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n34 <- intersect(C, D)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n134 <- intersect(n13, D)
      n234 <- intersect(n23, D)
      n1234 <- intersect(n123, D)
      grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
                                               area2 = length(B), area3 = length(C), area4 = length(D), 
                                               n12 = length(n12), n13 = length(n13), n14 = length(n14), 
                                               n23 = length(n23), n24 = length(n24), n34 = length(n34), 
                                               n123 = length(n123), n124 = length(n124), n134 = length(n134), 
                                               n234 = length(n234), n1234 = length(n1234), category = list.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, big.mark = ",",
                                               ...)
    }
    else if (5 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      E <- x[[5]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n15 <- intersect(A, E)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n25 <- intersect(B, E)
      n34 <- intersect(C, D)
      n35 <- intersect(C, E)
      n45 <- intersect(D, E)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n125 <- intersect(n12, E)
      n134 <- intersect(n13, D)
      n135 <- intersect(n13, E)
      n145 <- intersect(n14, E)
      n234 <- intersect(n23, D)
      n235 <- intersect(n23, E)
      n245 <- intersect(n24, E)
      n345 <- intersect(n34, E)
      n1234 <- intersect(n123, D)
      n1235 <- intersect(n123, E)
      n1245 <- intersect(n124, E)
      n1345 <- intersect(n134, E)
      n2345 <- intersect(n234, E)
      n12345 <- intersect(n1234, E)
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), 
                                                    area2 = length(B), area3 = length(C), area4 = length(D), 
                                                    area5 = length(E), n12 = length(n12), n13 = length(n13), 
                                                    n14 = length(n14), n15 = length(n15), n23 = length(n23), 
                                                    n24 = length(n24), n25 = length(n25), n34 = length(n34), 
                                                    n35 = length(n35), n45 = length(n45), n123 = length(n123), 
                                                    n124 = length(n124), n125 = length(n125), n134 = length(n134), 
                                                    n135 = length(n135), n145 = length(n145), n234 = length(n234), 
                                                    n235 = length(n235), n245 = length(n245), n345 = length(n345), 
                                                    n1234 = length(n1234), n1235 = length(n1235), 
                                                    n1245 = length(n1245), n1345 = length(n1345), 
                                                    n2345 = length(n2345), n12345 = length(n12345), 
                                                    category = list.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, ...)
      idx <- sapply(grob.list, function(i) grepl("text", i$name))
      
      for(i in 1:31){
        grob.list[idx][[i]]$label <- 
          format(as.numeric(grob.list[idx][[i]]$label), big.mark=",", scientific=FALSE)
      }
      
      grid.newpage()
      grid.draw(grob.list)
    }
    else {
      flog.error("Invalid size of input object", name = "VennDiagramLogger")
      stop("Invalid size of input object")
    }
  }
  if (length(x) == 2 & !is.null(total.population) & hyper.test) {
    val.p = calculate.overlap.and.pvalue(x[[1]], x[[2]], 
                                         total.population, lower.tail = lower.tail)
    if (is.null(sub)) {
      sub = paste0("p = ", signif(val.p[3], digits = 2))
    }
    else {
      sub = paste0(sub, ", p = ", signif(val.p[3], digits = 2))
    }
  }
  if (!is.null(sub)) {
    grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, 
                           fontface = sub.fontface, fontfamily = sub.fontfamily, 
                           col = sub.col, cex = sub.cex)
  }
  if (!is.null(main)) {
    grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, 
                           fontface = main.fontface, fontfamily = main.fontfamily, 
                           col = main.col, cex = main.cex)
  }
  if (!is.null(filename)) {
    current.type <- getOption("bitmapType")
    if (length(grep("Darwin", Sys.info()["sysname"]))) {
      options(bitmapType = "quartz")
    }
    else {
      options(bitmapType = "cairo")
    }
    if ("tiff" == imagetype) {
      tiff(filename = filename, height = height, width = width, 
           units = units, res = resolution, compression = compression)
    }
    else if ("png" == imagetype) {
      png(filename = filename, height = height, width = width, 
          units = units, res = resolution)
    }
    else if ("svg" == imagetype) {
      svg(filename = filename, height = height, width = width)
    }
    else {
      flog.error("You have misspelled your 'imagetype', please try again", 
                 name = "VennDiagramLogger")
      stop("You have misspelled your 'imagetype', please try again")
    }
    grid.draw(grob.list)
    dev.off()
    options(bitmapType = current.type)
    return(1)
  }
  return(grob.list)
}