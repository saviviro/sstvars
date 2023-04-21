#
wrap_multisets <- function(n, d) {
  # Laita kaksi accuracy vaihtoehtoa ja selvitä järkevä max dim?
  N <- choose(n + d - 1, d)
  if(N > 1e+8) {
    stop("Aborted since the large dimension required in the calculations may cause memory issues and likely take very long!")
  }
  get_multisets_Cpp(n=n, d=d, N=N)
}
