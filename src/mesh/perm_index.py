permutations = []
full_set = set(range(4))
for i in full_set:
  for j in full_set.difference({i}):
    for k in full_set.difference({i,j}):
      l = list(full_set.difference({i,j,k}))[0]
      permutations.append((i,j,k,l))

cmp_pairs = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]

def perm_idx(perm):
  n=0
  for a,b in cmp_pairs:
    n = 2 * n
    n += (perm[a] > perm[b])
  return n

i_perms = [(perm_idx(perm), perm) for perm in permutations]
i_perms.sort()
for p in i_perms:
  print(p)
    
