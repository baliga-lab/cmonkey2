/*
  gene_annotations.js - Javascript visualization library for MEME scores in
  a cMonkey cluster. This library is part of cmonkey2. See LICENSE for
  licensing details.
*/
var annot;
if (!annot) {
    annot = {};
}

(function () {
    "use strict";
    var pvalue_cutoff = 0.6;
    var base1_x = 500; // base 1 starts at this offset

    // draw a simple scale to see the sequence range size
    function drawScale(paper, params) {
        var tickLabels = ['-200', '-100', '-1'];
        var x = base1_x - 200, y = params.height - 30, width = 200, tickLen = 7, right = x + width, tickTop = y - tickLen, middle = x + width / 2;

        function drawTick(x, label) {
            var tick = paper.path('M' + x + ',' + tickTop + 'L' + x + ',' + y);
            var tickLabel = paper.text(x, y + 15, label);
        }
        var mainLine = paper.path('M' + x + ',' + y + 'L' + right + ',' + y);
        drawTick(x, tickLabels[0]);
        drawTick(middle, tickLabels[1]);
        drawTick(right, tickLabels[2]);
    }

    function drawCaption(paper, params) {
        var captions = ['log10(P) upstream meme', 'log10(P.clust)=NA; 8 seqs; 4 unique'];
        var x1 = 70, y1 = 10, x2 = 500, y2 = 10;
        var leftTop = paper.text(x1, y1, captions[0]);
        var rigthTop = paper.text(x2, y2, captions[1]);
    }

    function drawAnnotation(paper, annot, annotY, params) {
        var marginRight = 20, boxY = annotY, boxWidth = 100, boxHeight = 20, textOffsetY = boxHeight / 2,
        boxX = params.width - boxWidth - marginRight;
        var lineX = 40, lineY = boxY + boxHeight / 2;
        var annotHeight = 10;

        // box
        var box = paper.rect(boxX, boxY, boxWidth, boxHeight);
        box.attr('fill', annot.boxColor);
        box.attr('stroke', 'none');
        var geneCaption = paper.text(boxX + boxWidth / 2, boxY + textOffsetY, annot.gene);

        // annotation line
        var annotLine = paper.path('M' + lineX + ',' + lineY + 'L' + boxX + ',' + lineY);
        annotLine.attr('stroke', annot.lineColor);
        // p-value
        var pvalCaption = paper.text(lineX - 15, lineY, annot.log10.toString());

        // matches
        for (var i = 0; i < annot.matches.length; i++) {
            var match = annot.matches[i];
            //var matchX = lineX + match.start;
            var matchY = match.reverse ? lineY : lineY - annotHeight;
            var matchWidth = match.length;
            var matchX = base1_x - match.start - matchWidth;
            var matchbox = paper.rect(matchX, matchY, matchWidth, annotHeight);
            // p-values above the cutoff are transparent
            var opacity = match.score > pvalue_cutoff ? 0.0 : 1.0 - (match.score * 1.5);
            matchbox.attr('stroke', 'none');
            matchbox.attr('fill', match.motif == 0 ? '#f00' : '#0f0');
            matchbox.attr('fill-opacity', opacity.toString());
        }
        // gene description, optional
        //var condCaption = paper.text(lineX + 200, lineY, annot.condition);
    }

    function drawAnnotations(paper, params) {
        var annotations = params.annotations;
        var annotY = 20;
        for (var i = 0; i < annotations.length; i++) {
            drawAnnotation(paper, annotations[i], annotY, params);
            annotY += 25;
        }
    }

    function drawVerticals(paper, params) {
        var x1 = 400;
        var y1 = 0,  y2 = params.height;
        function drawVertical(x) {
            var vert = paper.path('M' + x + ',' + y1 + 'L' + x + ',' + y2);
            vert.attr('stroke-dasharray', '--');
            vert.attr('stroke', '#aaa');
        }
        //drawVertical(x1);
        drawVertical(base1_x);
    }

    annot.draw = function (id, params) {
        var paper = Raphael(document.getElementById(id),
                            params.width, params.height);
        var border = paper.rect(0, 0, params.width, params.height);
        base1_x = params.width - 140;
        drawVerticals(paper, params);
        drawScale(paper, params);
        //drawCaption(paper, params);
        drawAnnotations(paper, params);
    };
})();
