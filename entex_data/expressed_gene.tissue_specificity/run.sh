awk '$2==1' pg.txt | wc -l
awk '$2>=2' pg.txt | awk '$2<=5' | wc -l
awk '$2>=6' pg.txt | awk '$2<=25' | wc -l
awk '$2>25' pg.txt | wc -l
