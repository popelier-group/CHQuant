#Manually defined functional groups for molecules studied

molecule_atomcounts = {"alanine":22, "aspirin":21, 
              "azobenzene-cis":24, "azobenzene-trans":24,
              "glycine":19, "paracetamol":20, "urocanic-acid":16,"water":3}

molecule_functional_groups = {
"alanine":([[14,19,20,21],
                 [7,12,11,10],
                 [0,18,17,16],
                 [1,2],
                 [5,6]],
               ["Methyl", 
               "Methyl", 
               "Methyl", 
               "Carbonyl", 
               "Carbonyl"]),
    
"aspirin":([[12,17,18,19],
            [4,5,6,7,8,9],
            [2,10],
            [11,3]
           ],
           ["Methyl",
           "Ring",
           "Carbonyl",
           "Carbonyl"]
          ),
"azobenzene-cis":([[2,6,10,12,8,4],
                  [3,5,11,13,9,7],
                  [0,1]],
                 ["Ring",
                 "Ring",
                 "NN"]), #change "other"label to "Hs"; or add H's for each ring?
    
"azobenzene-trans":([[2,6,10,12,8,4],
                  [3,5,11,13,9,7],
                  [0,1]],
                 ["Ring",
                 "Ring",
                 "NN"]),

"glycine":([[11,16,18,17],
           [14,13,0,15],
           [5,6],
           [1,2],
           [4,8,9]#,
#            [10,12],
#            [3,7]
           ]
            ,
            ["Methyl",
            "Methyl",
            "Carbonyl",
            "Carbonyl",
            "Methylene"#,
#             "NH",
#             "NH"
            ]),

"paracetamol":([[10,17,16,18],
               [7,8,6,4,3,5],
               [1,9],
               [2,13],
               [12,15,11,14]#,
#                [0,19]
               ],
                 ["Methyl",
                  "Ring C",
                  "Carbonyl",
                  "NH",
                  "Ring H"#,
#                   "OH"
                 ]),

"urocanic-acid":([[7,3,6,4,2],
                 [13,10,12],
                 [0,15],
                 [1,9],
                 [5,8]],
                 ["Diazole C/N",
                 "Diazole H",
                 "OH",
                 "Carbonyl",
                 "C=C"]),

"water":([[0,1,2]],
         ["Water"]
        ),

}