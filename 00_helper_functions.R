`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
# function to create an interactive data table
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####----------------------------------------------------------------------#####
# convert a gene from all upper-case to the format in our data
#####----------------------------------------------------------------------#####

to_lower_gene_symbol <- function(s){
  list_string <- unlist(strsplit(s, ""))
  
  list_string[1] <- toupper(list_string[1])
  list_string[2:length(list_string)] <- tolower(list_string[2:length(list_string)])
  new.string <- paste(list_string, collapse = "")
  return(new.string)
}
