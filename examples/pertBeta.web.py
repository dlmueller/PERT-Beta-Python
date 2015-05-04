# -*- coding: utf-8 -*-
"""
Formulas from
    - Davis, R. 2008.
    - Teaching Note - Teaching Project Simulation in Excel Using PERT-Beta Distributions.
    - INFORMS Trans. Ed. 8(3) 139-148. Available online at http://ite.pubs.informs.org/.

@copyright: David Lukas Müller (2013, 2015)
"""

#--- infra
import os
import sys
import socket
import csv

#--- beta
import random
import math
import operator


#--- .
sys.path.append(os.path.join(os.getcwd(), '..'))
from pertbeta.betadist import BetaDistribution
from pertbeta.betadist import MAX_DOTS
from pertbeta import betadist

#--- webserver
import BaseHTTPServer

#---
def generateValues(dist, N) :
    values = [dist.random() for i in xrange(0, N)]
    return values

def makeBuckets(values, f) :
    buckets = {}
    valint = map(f, values)
    for v in valint :
        try :
            buckets[v] += 1
        except KeyError :
            buckets[v] = 1
    return buckets

def iterHistogram(buckets, maxDots = MAX_DOTS) :
    keys = sorted(buckets.keys())
    counts = buckets.values()
    N = sum(counts)
    bucketProps = dict([(v, buckets.get(v, 0) / float(N)) for v in keys])
    maxProp = max(bucketProps.values())
    scale = float(maxDots) / maxProp
    header = "%5s %6s %s" % ("units", "CDF", "PDF")
    yield header
    cdf = 0.0
    for x in keys :
        prop = bucketProps.get(x, 0.0)
        cdf += prop
        cdfPercent = 100.0 * cdf
        dots = round(scale * prop)
        theDots = "*"*dots
        line = "%(x)5i %(cdfPercent)6.2f %(theDots)s" % locals()
        yield line

def statValues(values) :
    count = float(len(values))
    mean = sum(values) / count
    var = sum([(v - mean) ** 2 for v in values]) / count
    return mean, var

def statValues2(values, confidence) :
    sortedValues = sorted(values)
    count = float(len(values))
    q25 = sortedValues[round(0.25 * count)]
    q50 = sortedValues[round(0.50 * count)]
    q75 = sortedValues[round(0.75 * count)]
    qConf = sortedValues[round(confidence * count)]
    return q25, q50, q75, qConf

def iter_properties(dist, a, m, b, alpha, beta) :
    yield "BETA(%(a).1f, %(m).1f, %(b).1f) alpha=%(alpha).1f beta=%(beta).1f" % locals()
    yield ""

def iter_analytisch(dist, a, m, b, alpha, beta) :
    mean = (a + 4 * m + b) / 6.0
    mean = a + (b - a) * (alpha / (alpha + beta))
    var = (b - a) ** 2 / 36.0
    var = (alpha / (alpha + beta)) * (beta / (alpha + beta)) * ((b - a) ** 2 / (alpha + beta + 1))
    sigma = math.sqrt(var)
    yield "ANALYTISCH"
    yield "mean1 = %(mean).1f var = %(var).1f sigma1 = %(sigma).1f" % locals()
    for line in dist.iterPDFasHistogram(epsilon = 0.0001) :
        yield line
    yield ""


def iter_empirisch(dist, a, m, b, alpha, beta, N) :
    if True :
        return # interessiert nicht mehr
    values = generateValues(dist, N)
    mean2, var2 = statValues(values)
    sigma2 = math.sqrt(var2)
    yield "EMPIRISCH"
    yield "mean2 = %(mean2).1f var = %(var2).1f sigma2 = %(sigma2).1f (N = %(N)i)" % locals()
    buckets = makeBuckets(values, round)
    for line in iterHistogram(buckets) :
        yield line
    yield ""

def iter_statBetaDist(a, m, b, N) :
    dist = BetaDistribution.FromAmB(a, m, b)
    alpha = dist.alpha
    beta = dist.beta
    for line in iter_properties(dist, a, m, b, alpha, beta):
        yield line
    for line in iter_analytisch(dist, a, m, b, alpha, beta) :
        yield line
    for line in iter_empirisch(dist, a, m, b, alpha, beta, N) :
        yield line

def statBetaDist(fout, a, m, b, N) :
    for line in iter_statBetaDist(a, m, b, N) :
        fout.write(line + "\n")
    return

def main(fout, param) :
    # interactive
    (o, n, p) = param
    fout.write("O/N/P = %(o)i/%(n)i/%(p)i<br/><br/>" % locals())
    N = 1000
    statBetaDist(fout, o, n, p, N)

def getParameterFromPath(thePath) :
    try :
        numbers = thePath.split("/")
        o = int(numbers[1])
        n = int(numbers[2])
        p = int(numbers[3])
        result = (o, n, p)
    except Exception as E :
        result = None
    return result

def getRandomParameters() :
    o = random.randint(1, 10)
    n = o + random.randint(1, 10)
    p = n + random.randint(1, 10)
    return (o, n, p)

def writeln(fout, line) :
    fout.write(line + "\n")

def writeHelp(fout, param):
    (o, n, p) = param
    writeln(fout, '<a href="http://localhost:8000/%(o)i/%(n)i/%(p)i">http://localhost:8000/%(o)i/%(n)i/%(p)i</a><br/><br/>' % locals())
    return

def iterBlogLink() :
    yield '<p>You can find lots of <a href="http://codingjourneyman.com/2014/10/06/the-clean-coder-estimation/">background information to PERT time estimations</a> in the web.</p>'
    yield '<p>PDF = Probability Density Function<br/>'
    yield 'CDF = Cumulative Distribution Function</p>'
    yield '<p>O = Optimistic Estimate<br/>'
    yield 'N = Nominal Estimate<br/>'
    yield 'P = Pessimistic Estimate</p>'

def iterUsage():
    yield '<pre>USAGE:<br>'
    yield 'python.exe pertBeta.web.py<br>'
    #yield 'python.exe pertBeta.web.py pertExample.csv<br>'
    yield '</pre>'
    for line in iterBlogLink() :
        yield line


def mathML_enabled():
    return False


def iterMathMlExample() :
    """
    Taken from http://www.xmlmind.com/tutorials/MathML/
    """
    if not mathML_enabled() :
        return
    yield "Beispiel f&uuml;r MathML"
    # Note the display="inline" attribute which specifies that the math element
    # is to be displayed inline (``like a word in a paragraph''). The other value
    # is display="block" which specifies that the math element is to be displayed
    # as a block (``like a paragraph''). This attribute has an influence on the
    # typographic rules applied by the MathML rendering engine.
    yield '<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline">'


    # mrow
    # Use this element to group any number of subexpressions horizontally.

    # mi
    # Use this element to specify an identifier, that is, the name of a variable,
    # a constant, a function, etc.
    # If this name is just one character long, the identifier is automatically
    # rendered using an italic font, otherwise the name is rendered using a normal,
    # upright, font.

    # mo
    # Use this element to specify an operator (e.g. '+'), a fence (e.g. '{') or
    # a separator (e.g. ',').
    # The appropriate amount of space is added on the left and on the right of an
    # mo depending on the textual contents of this element. Example: if in the
    # above expression you replace <mo>+</mo> by <mo>,</mo>, this will suppress
    # the space at the left of the mo element.

    # mn
    # Use this element to specify a numeric literal.
    yield '    <mfrac>'
    yield '  <mrow>'
    yield '    <msqrt>'
    yield '      <mi>x</mi>'
    yield '      <mo>+</mo>'
    yield '      <mi>y</mi>'
    yield '    </msqrt>'
    yield '  </mrow>'
    yield '    <mn>2</mn>'
    yield ' </mfrac>'
    yield '    <mo>=</mo>'
    yield '    <mn>2</mn>'
    yield '</math>'


def iterPageHeader():
    yield '<?xml version="1.0" encoding="UTF-8"?>'
    # next line will not work on INternet Explorer due to security considerations
    # http://www.matheboard.de/docs/W3C/Math/XSL/putting_math_on_the_web_german.htm
    # yield '<?xml-stylesheet type="text/xsl" href="http://www.w3.org/Math/XSL/mathml.xsl"?>'
    yield '<?xml-stylesheet type="text/xsl" href="mathml.xsl"?>'
    yield '<?xml-stylesheet type="text/xsl" href="pmathml.xsl"?>'
    # yield '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1 plus MathML 2.0//EN" "http://www.w3.org/Math/DTD/mathml2/xhtml-math11-f.dtd">'
    yield '<html xmlns="http://www.w3.org/1999/xhtml">'
    #yield '<html xmlns="http://www.w3.org/1999/xhtml" xmlns:xi="http://www.w3.org/2001/XInclude" xml:lang="de">'

def iterPageFooter() :
    yield '</html>'



def writeRequestContentLines(fout, theCommand, thePath, inputCsv) :
    """
    @param inputCsv: optionaler Name der Eingabedatei
    @type  inputCsv: C{str | None}
    """
    writeln(fout, "\n".join(iterPageHeader()))
    writeln(fout, "<body>")
    writeln(fout, "".join(iterUsage()))
    writeln(fout, "<pre>")
    writeln(fout, "%s %s" % (theCommand, thePath,))
    writeln(fout, "</pre>")

    param = getParameterFromPath(thePath)
    if param :
        write_REST_lines(fout, param, "REST")
    elif inputCsv :
        write_CSV_lines(fout, inputCsv, "CSV")
    else :
        param = getRandomParameters()
        write_REST_lines(fout, param, "Random")

    writeln(fout, "</body>")
    writeln(fout, "\n".join(iterPageFooter()))
    return

def write_CSV_lines(fout, inputCsv, heading) :
    writeln(fout, "<b>%(heading)s</b><br>" % locals())
    writeln(fout, "".join(iterInputCsv(inputCsv)))
    writeln(fout, "".join(iterMathMlExample()))

class MultiEstimate(object) :
    def __init__(self) :
        self._estimates = {} # ident -> dist

    def AppendEstimate(self, ident, dist):
        self._estimates[ident] = dist

    def getHeaderFields(self) :
        return ["ident",
                "opt", "likly", "pess",
                "alpha", "beta",
                "mean", "sigma",
                "mu+sig", "mu+2*sig",
                "curve",
                ]

    def getDataFields(self, ident, dist) :
        a = dist.a
        b = dist.b
        m = dist.m
        alpha = dist.alpha
        beta = dist.beta

        # Mittelwert und Streuung
        mean = dist.mean()
        sigma = dist.sigma()

        # Mittelwert plus Streuung
        mu_1sigma = mean + 1 * sigma
        mu_2sigma = mean + 2 * sigma

        # Kurve
        curve = '<img src="1pixel_green.png" />'

        return [ident,
                a, m, b,
                alpha, beta,
                mean, sigma,
                mu_1sigma, mu_2sigma,
                curve,
                ]

    def getFooterFields(self) :
        # Kennung
        ident = "ACCU"

        # O/N/P
        a = sum([dist.a for dist in self._estimates.values()])
        m = sum([dist.m for dist in self._estimates.values()])
        b = sum([dist.b for dist in self._estimates.values()])

        # alpha und beta
        alpha = ""
        beta = ""

        # Mittelwert und Streuung
        mean = sum([dist.mean() for dist in self._estimates.values()])
        sigma = math.sqrt(sum([dist.sigma() ** 2 for dist in self._estimates.values()]))

        # Mittelwert plus Streuung
        mu_1sigma = mean + 1 * sigma
        mu_2sigma = mean + 2 * sigma

        # Kurve
        curve = ""

        return [ident,
                a, m, b,
                alpha, beta,
                mean, sigma,
                mu_1sigma, mu_2sigma,
                curve,
                ]

    def cellFormat(self, pos) :
        """@rtype: C{str}"""
        colFormats = {1 : '.0f', # a
                      2 : '.0f', # m
                      3 : '.0f', # b
                      4 : '.3f', # alpha
                      5 : '.3f', # beta
                      6 : '.2f', # mean
                      7 : '.2f', # sigma
                      8 : '.1f', # mean + 1*sigma
                      9 : '.1f', # mean + 2*sigma
                      }
        theFormat = "%%%s" % (colFormats.get(pos, "s"),)
        return theFormat

    def cellHtmlAttributes(self, pos) :
        """@rtype: C{[str]}"""
        if pos in [1, 2, 3] :
            return ['align="center"']
        if pos in [4, 5, 6, 7] :
            return ['align="right"']
        if pos in [8, 9] :
            return ['align="right"']
        return []

    def formatCell(self, posContent) :
        """@rtype: C{str}"""
        (pos, content) = posContent
        if not content :
            return "&nbsp;"
        try :
            n = float(content)
            numberFormat = self.cellFormat(pos)
            return numberFormat % (n,)
        except ValueError as E:
            return content

    def iterFields(self):
        yield self.getHeaderFields()
        for (ident, dist) in self._estimates.iteritems() :
            yield map(self.formatCell, enumerate(self.getDataFields(ident, dist)))
        yield map(self.formatCell, enumerate(self.getFooterFields()))

def iterInputCsv(inputCsv) :
    with open(inputCsv, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ';', quotechar = '|')
        yield '<table border="1px">\n'
        me = MultiEstimate()
        for (rowNumber, rowList) in enumerate(spamreader):
            if rowNumber > 0 :
                a = int(rowList[0])
                m = int(rowList[1])
                b = int(rowList[2])
                ident = rowList[3]
                dist = BetaDistribution.FromAmB(a, m, b)
                me.AppendEstimate(ident, dist)
                continue

        for (row, colData) in enumerate(me.iterFields()) :
            columns = rowList + ["&nbsp;"]
            rowCells = ['<td %s>%s</td>' % (" ".join(me.cellHtmlAttributes(col)), s,)
                        for (col, s) in enumerate(map(str, colData))]
            row = "".join(rowCells)
            yield '<tr>%s</tr>\n' % (row,)
        yield '</table>\n'

def write_REST_lines(fout, param, heading):
    writeln(fout, "<b>%(heading)s</b><br>" % locals())
    writeln(fout, "<pre>")
    writeHelp(fout, param)
    main(fout, param)
    writeln(fout, "</pre>")
    writeln(fout, "".join(iterMathMlExample()))

def run(inputCsv):
    """
    @param inputCsv: optionaler Name der Eingabedatei
    @type  inputCsv: C{str | None}
    """

    class RequestHandler(BaseHTTPServer.BaseHTTPRequestHandler) :

        def do_GET(self) :
            # if self.path.endswith('.html'):
            # content = "<html><body>hallo</body></html>\n"
            self.send_response(200)
            self.send_header('Content-type', 'text/html')
            # self.send_header('Content-Length', len(content))
            self.end_headers()
            writeRequestContentLines(self.wfile, self.command, self.path, inputCsv)
            return

    port = 8000
    targetUrl = "http://localhost:%(port)i/" % locals()
    server_address = ('', port)
    httpd = BaseHTTPServer.HTTPServer(server_address, RequestHandler)
    os.startfile(targetUrl)
    httpd.serve_forever()

if __name__ == "__main__" :
    try :
        inputCsv = sys.argv[1]
    except IndexError :
        inputCsv = None
    run(inputCsv)
