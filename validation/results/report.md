| method | stratum | queries | accepted | true | false | false-assoc. rate (95% CI) | FP/1,000 | sensitivity |
|---|---|---|---|---|---|---|---|---|
| method1_naive | all | 1000 | 909 | 179 | 730 | 80.3% (77.6-82.8) | 730.0 |  |
| method1_naive | has_structure | 215 | 213 | 176 | 37 | 17.4% (12.9-23.0) | 172.1 | 81.9% |
| method1_naive | no_structure | 785 | 696 | 3 | 693 | 99.6% (98.7-99.9) | 882.8 |  |
| method1_validated_strict | all | 1000 | 85 | 85 | 0 | 0.0% (0.0-4.3) | 0.0 |  |
| method1_validated_strict | has_structure | 215 | 85 | 85 | 0 | 0.0% (0.0-4.3) | 0.0 | 39.5% |
| method1_validated_strict | no_structure | 785 | 0 | 0 | 0 | n/a | 0.0 |  |
| method1_validated_permissive | all | 1000 | 160 | 105 | 55 | 34.4% (27.5-42.0) | 55.0 |  |
| method1_validated_permissive | has_structure | 215 | 109 | 104 | 5 | 4.6% (2.0-10.3) | 23.3 | 48.4% |
| method1_validated_permissive | no_structure | 785 | 51 | 1 | 50 | 98.0% (89.7-99.7) | 63.7 |  |
| method2_naive | all | 1000 | 343 | 189 | 154 | 44.9% (39.7-50.2) | 154.0 |  |
| method2_naive | has_structure | 215 | 210 | 186 | 24 | 11.4% (7.8-16.4) | 111.6 | 86.5% |
| method2_naive | no_structure | 785 | 133 | 3 | 130 | 97.7% (93.6-99.2) | 165.6 |  |
| method2_validated_strict | all | 1000 | 85 | 85 | 0 | 0.0% (0.0-4.3) | 0.0 |  |
| method2_validated_strict | has_structure | 215 | 85 | 85 | 0 | 0.0% (0.0-4.3) | 0.0 | 39.5% |
| method2_validated_strict | no_structure | 785 | 0 | 0 | 0 | n/a | 0.0 |  |
| method2_validated_permissive | all | 1000 | 160 | 105 | 55 | 34.4% (27.5-42.0) | 55.0 |  |
| method2_validated_permissive | has_structure | 215 | 108 | 104 | 4 | 3.7% (1.5-9.1) | 18.6 | 48.4% |
| method2_validated_permissive | no_structure | 785 | 52 | 1 | 51 | 98.1% (89.9-99.7) | 65.0 |  |
| method3_strict | all | 1000 | 83 | 83 | 0 | 0.0% (0.0-4.4) | 0.0 |  |
| method3_strict | has_structure | 215 | 83 | 83 | 0 | 0.0% (0.0-4.4) | 0.0 | 38.6% |
| method3_strict | no_structure | 785 | 0 | 0 | 0 | n/a | 0.0 |  |
| method3_permissive | all | 1000 | 144 | 96 | 48 | 33.3% (26.2-41.4) | 48.0 |  |
| method3_permissive | has_structure | 215 | 97 | 95 | 2 | 2.1% (0.6-7.2) | 9.3 | 44.2% |
| method3_permissive | no_structure | 785 | 47 | 1 | 46 | 97.9% (88.9-99.6) | 58.6 |  |

False associations by tier (all queries):

| method | same sequence | near-identical | different protein | unmeasured |
|---|---|---|---|---|
| method1_naive | 34 | 70 | 503 | 123 |
| method1_validated_strict | 0 | 0 | 0 | 0 |
| method1_validated_permissive | 34 | 21 | 0 | 0 |
| method2_naive | 34 | 75 | 15 | 30 |
| method2_validated_strict | 0 | 0 | 0 | 0 |
| method2_validated_permissive | 34 | 21 | 0 | 0 |
| method3_strict | 0 | 0 | 0 | 0 |
| method3_permissive | 31 | 17 | 0 | 0 |

Correct naive method-1 hits, by what the strict gate did with them:

| outcome | n |
|---|---|
| accepted (M4) | 85 |
| partial structure (accession confirmed, identity >= 95, coverage < 80) | 35 |
| unmeasured (no coordinate file validated) | 23 |
| M3 (high identity and coverage, accession not confirmed) | 20 |
| other (M1) | 12 |
| other (M2) | 4 |
