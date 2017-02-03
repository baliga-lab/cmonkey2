/* seqlogo.js - see README and LICENSE for details */
var seqlogo;
if (!seqlogo) {
    seqlogo = {};
}
(function () {
    "use strict";
    // some default settings
    var MARGIN_LEFT = 40, MARGIN_TOP = 20, MARGIN_RIGHT = 20,
        MARGIN_BOTTOM = 30, DEFAULT_OPTIONS, NUCLEOTIDE_COLORS,
        AMINO_COLORS, MEASURE_CANVAS, STRETCH = 0.65, BASELINE = 6;
    NUCLEOTIDE_COLORS = {
        'A': 'rgb(0, 128, 0)',
        'G': 'rgb(255, 165, 0)',
        'T': 'rgb(255, 0, 0)',
        'U': 'rgb(255, 0, 0)',
        'C': 'rgb(0, 0, 255)'
    };
    AMINO_COLORS = {
        // polar amino acids
        'G': 'rgb(0, 200, 50)',
        'S': 'rgb(0, 200, 50)',
        'T': 'rgb(0, 200, 50)',
        'Y': 'rgb(0, 200, 50)',
        'C': 'rgb(0, 200, 50)',
        'Q': 'rgb(0, 200, 50)',
        'N': 'rgb(0, 200, 50)',
        // basic
        'K': 'rgb(0, 0, 230)',
        'R': 'rgb(0, 0, 230)',
        'H': 'rgb(0, 0, 230)',
        // acidic
        'D': 'rgb(255, 0, 0)',
        'E': 'rgb(255, 0, 0)',
        // hydrophobic
        'A': 'rgb(0, 0, 0)',
        'V': 'rgb(0, 0, 0)',
        'L': 'rgb(0, 0, 0)',
        'I': 'rgb(0, 0, 0)',
        'P': 'rgb(0, 0, 0)',
        'W': 'rgb(0, 0, 0)',
        'F': 'rgb(0, 0, 0)',
        'M': 'rgb(0, 0, 0)'
    };
    DEFAULT_OPTIONS = {
        type: 'canvas',
        width: 400,
        height: 300,
        glyphStyle: '20pt Helvetica'
    };
    MEASURE_CANVAS = document.createElement('canvas');
    MEASURE_CANVAS.setAttribute('width', 500);
    MEASURE_CANVAS.setAttribute('height', 500);

    // **********************************************************************
    // ****** Common Functions
    // **********************************************************************

    function rank(arr) {
        var result = [], i;
        for (i = 0; i < arr.length; i += 1) {
            result.push([i, arr[i]]);
        }
        return result.sort(function (a, b) {
            return a[1] - b[1];
        });
    }

    function log(n, base) {
        return Math.log(n) / Math.log(base);
    }
    function uncertaintyAt(pssm, motifPos) {
        var row, freq, sum = 0;
        for (row = 0; row < pssm.values[motifPos].length; row += 1) {
            freq = pssm.values[motifPos][row];
            if (freq > 0) {
                sum +=  freq * log(freq, 2);
            }
        }
        return -sum;
    }
    function rsequence(pssm, motifPos) {
        var correctionFactor = 0.0, numBits;
        numBits = Math.ceil(log(pssm.alphabet.length, 2));
        return numBits - (uncertaintyAt(pssm, motifPos) + correctionFactor);
    }

    // Generic PSSM drawing function
    function drawPSSM(pssm, scalex, y0, yHeight, drawFun) {
        var x, y, motifPos, size, columnRanks, currentGlyph, row, maxWidth, rseq;
        x = MARGIN_LEFT;

        for (motifPos = 0; motifPos < pssm.values.length; motifPos += 1) {
            y = y0;
            columnRanks = rank(pssm.values[motifPos]);
            maxWidth = 0;
            rseq = rsequence(pssm, motifPos);
            for (row = 0; row < columnRanks.length; row += 1) {
                currentGlyph = pssm.alphabet[columnRanks[row][0]];
                size = drawFun(currentGlyph, x, y, scalex, yHeight,
                               rseq * columnRanks[row][1]);
                if (size.width > maxWidth) {
                    maxWidth = size.width;
                }
                y -= size.height;
            }
            x += maxWidth;
        }
        return x;
    }

    // **********************************************************************
    // ****** Canvas-based Implementation
    // **********************************************************************

    function firstLine(imageData) {
        var pixels = imageData.data, row, col, index;
        for (row = 0; row < imageData.height; row += 1) {
            for (col = 0; col < imageData.width; col += 1) {
                index = (row * imageData.width * 4) + col * 4;
                if (pixels[index + 3] !== 0) {
                    return row;
                }
            }
        }
        return imageData.height;
    }
    function lastLine(imageData) {
        var pixels = imageData.data, row, col, index;
        for (row = imageData.height - 1; row >= 0; row -= 1) {
            for (col = 0; col < imageData.width; col += 1) {
                index = (row * imageData.width * 4) + col * 4;
                if (pixels[index + 3] !== 0) {
                    return row;
                }
            }
        }
        return imageData.height - 1;
    }

    function measureText(text, font, scalex, scaley) {
        var imageData, context, first;
        if (scaley === 0) {
            return 0;
        }
        context = MEASURE_CANVAS.getContext('2d');
        context.fillStyle = "rgb(0, 0, 0)";
        context.font = font;
        context.textBaseline = 'top';
        context.save();
        context.scale(scalex, scaley);
        context.fillText(text, 0, 0);
        context.restore();

        imageData = context.getImageData(0, 0, MEASURE_CANVAS.width, MEASURE_CANVAS.height);
        first = firstLine(imageData);
        context.clearRect(0, 0, MEASURE_CANVAS.width, MEASURE_CANVAS.height);
        return lastLine(imageData) - first + 1;
    }

    function drawLabelsY(context, numBits, x0, y0, yHeight) {
        var i, label, ydist = yHeight / numBits, y = y0;

        context.font = '12pt Arial';
        context.fillText('bits', x0 + 10, MARGIN_TOP - 5);
        var textHeight = measureText('M', context.font, 1.0, 1.0);
        y += textHeight / 2;

        for (i = 0; i <= numBits; i += 1) {
            label = i.toString();
            context.fillText(label, x0, y);
            y -= ydist;
        }
    }

    function drawMinorTicksY(context, y0, y1, numDivisions) {
        var interval = (y1 - y0) / numDivisions, y = y0;
        for (var i = 0; i < numDivisions; i++) {
            if (i > 0) {
                context.beginPath();
                context.moveTo(MARGIN_LEFT - 5, y);
                context.lineTo(MARGIN_LEFT, y);
                context.stroke();
            }
            y += interval;
        }
    }

    function drawTicksY(context, numBits, bottom) {
        var mainIntervalY = (bottom - MARGIN_TOP) / numBits;
        var y = MARGIN_TOP;
        for (var i = 0; i <= numBits; i++) {
            context.beginPath();
            context.moveTo(MARGIN_LEFT - 10, y);
            context.lineTo(MARGIN_LEFT, y);
            context.stroke();
            if (i < numBits) drawMinorTicksY(context, y, y + mainIntervalY, 5);
            y += mainIntervalY;
        }
    }

    function drawAxis(context, numBits, right, bottom) {
        // main axis
        context.beginPath();
        context.moveTo(MARGIN_LEFT, MARGIN_TOP);
        context.lineTo(MARGIN_LEFT, bottom);
        context.lineTo(right, bottom);
        context.stroke();

        drawTicksY(context, numBits, bottom);
    }

    function drawScale(canvas, pssm) {
        var context = canvas.getContext('2d'), right = canvas.width - MARGIN_RIGHT,
            numBits = Math.ceil(log(pssm.alphabet.length, 2)),
            bottom = canvas.height - MARGIN_BOTTOM;
        drawAxis(context, numBits, right, bottom);
        drawLabelsY(context, numBits, MARGIN_LEFT - 25, bottom, bottom - MARGIN_TOP);
    }

    function drawGlyph(context, glyph, colors, x, y, scalex,
                       yHeight, maxFontHeightNormal,
                       weight) {
        var glyphWidth, scaley, glyphHeightScaled;
        glyphWidth = context.measureText(glyph).width * scalex;
        scaley = weight * (yHeight / maxFontHeightNormal) * STRETCH;
        glyphHeightScaled = measureText(glyph, context.font, scalex, scaley);
        if (scaley > 0) {
            context.fillStyle = colors[glyph];
            context.save();
            context.translate(x, y);
            context.scale(scalex, scaley);
            context.translate(-x, -y);
            context.fillText(glyph, x, y);
            context.restore();
        }
        return { width: glyphWidth, height: glyphHeightScaled };
    }

    function colorTableFor(pssm) {
        var i, c;
        for (i = 0; i < pssm.alphabet.length; i += 1) {
            c = pssm.alphabet[i];
            if (c !== 'A' && c !== 'G' && c !== 'T' && c !== 'C' && c !== 'U') {
                return AMINO_COLORS;
            }
        }
        return NUCLEOTIDE_COLORS;
    }

    function drawGlyphs(canvas, options, pssm) {
        var context, yHeight, maxFontHeightNormal, sumColumnWidthsNormal, xWidth, scalex;
        context = canvas.getContext('2d');
        context.textBaseline = 'alphabetic';
        context.font = options.glyphStyle;
        yHeight = canvas.height - (MARGIN_BOTTOM + MARGIN_TOP);
        maxFontHeightNormal = measureText('Mg', context.font, 1.0, 1.0);
        sumColumnWidthsNormal = context.measureText('W').width * pssm.values.length;
        xWidth = canvas.width - (MARGIN_LEFT + MARGIN_RIGHT);
        scalex = xWidth / sumColumnWidthsNormal;
        var lastX = drawPSSM(pssm, scalex,
                             canvas.height - MARGIN_BOTTOM, yHeight,
                             function (currentGlyph, x, y, scalex, yHeight, weight) {
                                 return drawGlyph(context, currentGlyph,
                                                  colorTableFor(pssm), x, y,
                                                  scalex, yHeight,
                                                  maxFontHeightNormal, weight);
                             });
        return (lastX - MARGIN_LEFT) / pssm.values.length;
    }

    function drawTicksX(canvas, pssm,  interval) {
        var context = canvas.getContext('2d'), bottom = canvas.height - MARGIN_BOTTOM;
        context.font = '12pt Arial';
        context.fillStyle = 'black';
        for (var i = 1; i <= pssm.values.length; i++) {
            var x = MARGIN_LEFT + i * interval;
            var xi = x - interval / 2;
            var tickHeight = (i % 5 == 0) ? 10 : 5;
            context.beginPath();
            context.moveTo(xi, bottom);
            context.lineTo(xi, bottom + tickHeight);
            context.stroke();
            if (i % 5 == 0) {
                var label = i.toString();
                var textdim = context.measureText(label);
                // the TextMetrics object currently does not have any other attributes
                // than width, so we simply specify a text height
                var textHeight = 14;
                context.fillText(label, xi - textdim.width / 2, bottom + tickHeight + textHeight);
            }
        }
    }

    function makeCanvas(id, options, pssm) {
        var canvas = document.createElement("canvas"), elem;
        canvas.id = id;
        canvas.setAttribute('width', options.width);
        canvas.setAttribute('height', options.height);
        canvas.setAttribute('style', 'border: 1px solid black');
        elem = document.getElementById(id);
        elem.parentNode.replaceChild(canvas, elem);
        drawScale(canvas, pssm);
        var interval = drawGlyphs(canvas, options, pssm);
        drawTicksX(canvas, pssm, interval);
    }

    // **********************************************************************
    // ****** Public API
    // **********************************************************************

    seqlogo.makeLogo = function (id, pssm, options) {
        if (options === null) {
            options = DEFAULT_OPTIONS;
        }
        // TODO: copy the options from DEFAULT_OPTIONS that are missing
        makeCanvas(id, options, pssm);
    };
}());
