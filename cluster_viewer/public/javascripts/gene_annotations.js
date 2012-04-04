var annot;
if (!annot) {
    annot = {};
}
(function () {
    "use strict";

    // draw a simple scale to see the sequence range size
    function drawScale(paper, params) {
        var tickLabels = ['-200', '-100', '-1'];
        var x = 300, y = params.height - 30, width = 200, tickLen = 7, right = x + width, tickTop = y - tickLen, middle = x + width / 2;

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
        var annotHeight = 15;

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
            var matchX = lineX + match.start, matchY = match.reverse ? lineY : lineY - annotHeight;
            var matchWidth = match.length;
            var matchbox = paper.rect(matchX, matchY, matchWidth, annotHeight);
            matchbox.attr('stroke', 'none');
            matchbox.attr('fill', match.motif == 0 ? '#f00' : '#0f0');
            matchbox.attr('fill-opacity', match.score.toString());
        }
        // gene description, optional
        //var condCaption = paper.text(lineX + 200, lineY, annot.condition);
    }

    function drawAnnotations(paper, params) {
        var annotations = params.annotations;
        var annotY = 50;
        for (var i = 0; i < annotations.length; i++) {
            drawAnnotation(paper, annotations[i], annotY, params);
            annotY += 40;
        }
    }

    function drawVerticals(paper, params) {
        var x1 = 400, x2 = 420;
        var y1 = 0,  y2 = params.height;
        function drawVertical(x) {
            var vert = paper.path('M' + x + ',' + y1 + 'L' + x + ',' + y2);
            vert.attr('stroke-dasharray', '--');
            vert.attr('stroke', '#aaa');
        }
        drawVertical(x1);
        drawVertical(x2);
    }

    annot.draw = function (id, params) {
        var paper = Raphael(document.getElementById(id),
                            params.width, params.height);
        var border = paper.rect(0, 0, params.width, params.height);
        drawVerticals(paper, params);
        drawScale(paper, params);
        //drawCaption(paper, params);
        drawAnnotations(paper, params);
    };
})();
