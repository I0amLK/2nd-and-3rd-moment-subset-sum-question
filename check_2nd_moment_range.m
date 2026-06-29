// Helper function to check the theorem for a single q and fixed k.
// Returns true and an empty string if all pairs are covered,
// otherwise false and an error message.
CheckTheoremForQ := function(q, k)
    F := GF(q);
    
    if Characteristic(F) lt 5 then
        return false, "Characteristic < 5";
    end if;
    if q lt 17 then
        return false, "q < 17";
    end if;
    if k lt 4 or k gt q-4 then
        return false, "k not in [4, q-4]";
    end if;
    
    S := {x : x in F};
    all_pairs := {<u,v> : u in F, v in F};
    covered := {};
    
    total_subsets := Binomial(q, k);
    printf "  q=%o, k=%o → %o subsets. ", q, k, total_subsets;
    
    if total_subsets gt 500000 then
        printf "WARNING: Very large search, may hang.\n";
        // Uncomment the next line to skip dangerous cases:
        // return false, "Too many subsets";
    end if;
    
    for I in Subsets(S, k) do
        sum1 := &+I;
        sum_sq := &+[x^2 : x in I];
        sum2 := (sum1^2 - sum_sq) / 2;
        Include(~covered, <sum1, sum2>);
        if #covered eq #all_pairs then
            printf "All %o pairs covered.\n", #all_pairs;
            return true, "";
        end if;
    end for;
    
    missing := all_pairs diff covered;
    miss_seq := SetToSequence(missing);
    first_few := miss_seq[1..Minimum(3, #miss_seq)];
    printf "Missing %o pairs, e.g. %o\n", #missing, first_few;
    return false, "";
end function;

// ============================================
//  Set the parameters for your run
// ============================================
k := 6;         // fixed subset size (4 <= k <= q-4)
Qmax := 31;     // maximal field size to check

printf "Checking theorem for k=%o and all admissible q <= %o\n", k, Qmax;
printf "(q must be a prime power, char >=5, q>=17, and 4<=k<=q-4)\n\n";

all_ok := true;
for q := 17 to Qmax do
    // only prime powers with suitable characteristic
    if not IsPrimePower(q) then continue; end if;
    fact := Factorization(q);
    p := fact[1,1];
    if p lt 5 then continue; end if;
    
    // check range condition
    if k lt 4 or k gt q-4 then
        printf "q=%o: k out of allowed range, skipping.\n", q;
        continue;
    end if;
    
    printf "Testing q=%o...\n", q;
    ok, err := CheckTheoremForQ(q, k);
    if not ok then
        printf "FAILED for q=%o: %o\n\n", q, err;
        all_ok := false;
        break;   // stop on first failure
    end if;
end for;

if all_ok then
    printf "\nAll tested q passed the theorem.\n";
end if;