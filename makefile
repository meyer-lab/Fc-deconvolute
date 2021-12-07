flist = 1 3 7 11
flistPath = $(patsubst %, output/figure%.svg, $(flist))

all: $(flistPath)

venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install --prefer-binary -Uqr requirements.txt
	touch venv/bin/activate

output/figure%.svg: venv genFigure.py deconv/figures/figure%.py
	@mkdir -p output
	. venv/bin/activate && ./genFigure.py $*

test: venv
	. venv/bin/activate && pytest -s -v -x

clean:
	rm -rf output venv
