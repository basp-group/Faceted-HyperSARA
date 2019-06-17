function id = interleaved_facets(L, Q)

if Q > L
    error('Number of facets Q=%i greater than the dimension L=%i', Q, L);
end

id = cell(Q, 1);
for q = 1:Q
    id{q} = q:Q:L;
end

end