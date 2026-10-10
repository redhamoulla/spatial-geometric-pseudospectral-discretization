"""Static, shareable SVG figures using ReportLab's scientific line-plot objects."""
import math
from reportlab.graphics.shapes import Drawing, String, Rect, Line
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics import renderSVG
from reportlab.lib.colors import HexColor, white

COLORS=["#087f8c","#c44c35","#5763ac","#bd8b24"]


def plot(path,title,subtitle,series,xlabel,ylabel):
    drawing=Drawing(780,450)
    drawing.add(Rect(0,0,780,450,fillColor=white,strokeColor=None))
    drawing.add(String(28,416,title,fontName="Helvetica-Bold",fontSize=20,fillColor=HexColor("#172c41")))
    drawing.add(String(28,392,subtitle,fontName="Helvetica",fontSize=10,fillColor=HexColor("#526475")))
    chart=LinePlot()
    chart.x,chart.y,chart.width,chart.height=80,112,660,255
    chart.data=[[(float(x),float(y)) for x,y in points if math.isfinite(x) and math.isfinite(y)] for _,points in series]
    chart.joinedLines=1
    chart.xValueAxis.labels.fontSize=9
    chart.yValueAxis.labels.fontSize=9
    chart.xValueAxis.strokeColor=chart.yValueAxis.strokeColor=HexColor("#8b9daa")
    chart.yValueAxis.visibleGrid=1
    chart.yValueAxis.gridStrokeColor=HexColor("#e2e9ee")
    chart.yValueAxis.gridStrokeWidth=.5
    chart.xValueAxis.labelTextFormat=lambda x:f"{x:.3g}"
    chart.yValueAxis.labelTextFormat=lambda x:f"{x:.6g}"
    for j in range(len(series)):
        chart.lines[j].strokeColor=HexColor(COLORS[j%len(COLORS)])
        chart.lines[j].strokeWidth=1.8
    drawing.add(chart)
    drawing.add(String(410,80,xlabel,textAnchor="middle",fontName="Helvetica",fontSize=11))
    drawing.add(String(30,378,ylabel,fontName="Helvetica",fontSize=10))
    for j,(name,_) in enumerate(series):
        xx=40+(j%2)*375; yy=53-(j//2)*20
        drawing.add(Line(xx,yy,xx+22,yy,strokeColor=HexColor(COLORS[j%4]),strokeWidth=2.5))
        drawing.add(String(xx+30,yy-3,name,fontName="Helvetica",fontSize=10))
    renderSVG.drawToFile(drawing,str(path))
    return drawing
