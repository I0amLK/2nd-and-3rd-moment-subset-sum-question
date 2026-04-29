CheckTheorem := function(q, k)
    F := GF(q);
    
    
    S := {x : x in F};
    all_pairs := {<u,v> : u in F, v in F};
    covered := {};
    
    printf "Checking all %o-subsets of F_%o (total %o subsets)\n", k, q, Binomial(q, k);
    
    for I in Subsets(S, k) do
        sum1 := &+I;
        sum_sq := &+[x^2 : x in I];
        // Division by 2 in the finite field (char ≠ 2)
        sum2 := (sum1^2 - sum_sq) / 2;
        Include(~covered, <sum1, sum2>);
        
        if #covered eq #all_pairs then
            printf "All %o pairs covered.\n", #all_pairs;
            return true;
        end if;
    end for;
    
    missing := all_pairs diff covered;
    if #missing eq 0 then
        return true;
    else
        miss_seq := SetToSequence(missing);
        first_few := miss_seq[1..Minimum(10, #miss_seq)];
        printf "Missing %o pairs. First few: %o\n", #missing, first_few;
        return false;
    end if;
end function;

CheckTheorem(17,4);