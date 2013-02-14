all: asms.txt \
	word-count.txt

COMPONENTS = title.txt \
	novel.txt \
	intro.txt \
	method.txt \
	data.txt


asms.txt: $(COMPONENTS) \
	word-count.txt
	echo "TITLE [max 20 words]: " > asms.txt
	cat title.txt >> asms.txt
	echo "" >> asms.txt
	echo "=======================================================" >> asms.txt
	echo "" >> asms.txt
	echo "AUTHORS: " >> asms.txt
	echo "Benjamin Diament, William Stafford Noble" >> asms.txt
	echo "" >> asms.txt
	echo "=======================================================" >> asms.txt
	echo "" >> asms.txt
	echo "NOVEL ASPECT [max 20 words]: " >> asms.txt
	cat novel.txt >> asms.txt
	echo "" >> asms.txt
	echo "=======================================================" >> asms.txt
	echo "" >> asms.txt
	echo "INTRODUCTION [max 120 words]: " >> asms.txt
	cat intro.txt >> asms.txt
	echo "" >> asms.txt
	echo "=======================================================" >> asms.txt
	echo "" >> asms.txt
	echo "METHOD [max 120 words]: " >> asms.txt
	cat method.txt >> asms.txt
	echo "" >> asms.txt
	echo "=======================================================" >> asms.txt
	echo "" >> asms.txt
	echo "PRELIMINARY DATA [max 300 words]: " >> asms.txt
	cat data.txt >> asms.txt
	echo "" >> asms.txt
	echo "=======================================================" >> asms.txt
	echo "" >> asms.txt
	cat word-count.txt >> asms.txt

word-count.txt: $(COMPONENTS)
	wc -w $(COMPONENTS) > word-count.txt

tar: $(COMPONENTS) makefile asms.txt word-count.txt
	tar cvf asms-2009.tar $^

clean:
	rm -f asms.txt word-count.txt asms-2009.tar
