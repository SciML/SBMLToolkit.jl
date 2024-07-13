a = first(@species A [isconstantspecies = false])
b = first(@species A [isbcspecies = true])
isequal(a, b)
