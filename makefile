# Build any figure files you find
flist = $(wildcard deconv/figures/figure*.py)
all: $(patsubst deconv/figures/figure%.py, output/figure%.svg, $(flist))

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
