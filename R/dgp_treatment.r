#-------------------------------------------------------------------
# %#     Random Treatment Assignment
#-------------------------------------------------------------------
gen.treat.sreg <- function(pi.matr.w, ns, k)
#-------------------------------------------------------------------
{
  rows <- nrow(pi.matr.w)
  code.elements <- character(rows + 1)

  for (i in 1:rows)
  {
    code.elements[i] <- paste0(
      "rep(", i,
      ", floor(pi.matr.w[", i, ",", k, "]*ns))"
    )
  }

  code.elements[rows + 1] <- paste0(
    "rep(0, (ns - ",
    paste0("floor(pi.matr.w[", 1:rows,
      ",", k, "]*ns)",
      collapse = " - "
    ), "))"
  )

  code <- paste(code.elements, collapse = ", ")

  result <- eval(parse(text = paste("sample(c(", code, "))")))

  return(result)
}
#-------------------------------------------------------------------
gen.treat.creg <- function(pi.matr.w, ns, k)
#-------------------------------------------------------------------
{
  rows <- nrow(pi.matr.w)
  code.elements <- character(rows + 1)

  for (i in 1:rows)
  {
    code.elements[i] <- paste0(
      "rep(", i,
      ", floor(pi.matr.w[", i, ",", k, "]*ns))"
    )
  }

  code.elements[rows + 1] <- paste0(
    "rep(0, (ns - ",
    paste0("floor(pi.matr.w[", 1:rows,
      ",", k, "]*ns)",
      collapse = " - "
    ), "))"
  )

  code <- paste(code.elements, collapse = ", ")

  result <- eval(parse(text = paste("sample(c(", code, "))")))

  return(result)
}
