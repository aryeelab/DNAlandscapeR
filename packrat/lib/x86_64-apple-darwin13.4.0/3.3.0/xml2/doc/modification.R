## ------------------------------------------------------------------------
library(xml2)
x <- read_xml("<p>This is some <b>text</b>. This is more.</p>")
xml_text(x)

xml_text(x) <- "This is some other text."
xml_text(x)

# You can avoid this by explicitly selecting the text node.
x <- read_xml("<p>This is some text. This is <b>bold!</b></p>")
text_only <- xml_find_all(x, "//text()")

xml_text(text_only) <- c("This is some other text. ", "Still bold!")
xml_text(x)
xml_structure(x)

## ------------------------------------------------------------------------
x <- read_xml("<a href='invalid!'>xml2</a>")
xml_attr(x, "href")

xml_attr(x, "href") <- "https://github.com/hadley/xml2"
xml_attr(x, "href")

xml_attrs(x) <- c(id = "xml2", href = "https://github.com/hadley/xml2")
xml_attrs(x)
cat(as.character(x))

xml_attrs(x) <- NULL
cat(as.character(x))

## ------------------------------------------------------------------------
x <- read_xml("<a><b/></a>")
x
xml_name(x)
xml_name(x) <- "c"
x

## ------------------------------------------------------------------------
x <- read_xml("<parent><child>1</child><child>2<child>3</child></child></parent>")
children <- xml_children(x)
t1 <- children[[1]]
t2 <- children[[2]]
t3 <- xml_children(children[[2]])[[1]]

xml_replace(t1, t3)
x

## ------------------------------------------------------------------------
x <- read_xml("<parent><child>1</child><child>2<child>3</child></child></parent>")
children <- xml_children(x)
t1 <- children[[1]]
t2 <- children[[2]]
t3 <- xml_children(children[[2]])[[1]]

xml_add_sibling(t1, t3)
x

xml_add_sibling(t3, t1, where = "before")
x

## ------------------------------------------------------------------------
x <- read_xml("<parent><child>1</child><child>2<child>3</child></child></parent>")
children <- xml_children(x)
t1 <- children[[1]]
t2 <- children[[2]]
t3 <- xml_children(children[[2]])[[1]]

xml_add_child(t1, t3)
x

xml_add_child(t1, read_xml("<test/>"))
x

## ------------------------------------------------------------------------
library(magrittr)
d <- xml_new_document() %>%
  xml_add_child("sld",
    xmlns = "http://www.o.net/sld",
    "xmlns:ogc" = "http://www.o.net/ogc",
    "xmlns:se" = "http://www.o.net/se",
    version = "1.1.0") %>%
  xml_add_child("layer") %>%
  xml_add_child("se:Name") %>%
  xml_root()

cat(as.character(d))

