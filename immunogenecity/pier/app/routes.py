from flask import render_template,flash,redirect,url_for,request
from app import app
from app.program import computing_s,computing_m,hla_convert,svg_path
from app.form import QueryForm

@app.route('/result')
def result():
    peptide = request.args.get('peptide')
    mhc = request.args.get('mhc')
    score = computing_s(peptide,mhc)
    p,m,i = computing_m(peptide,mhc)
    m3 = hla_convert(mhc)
    path = ["/static/{0}_positive_9.png".format(m3),"/static/{0}_negative_9.png".format(m3),"/static/{0}_positive_10.png".format(m3),"/static/{0}_negative_9.png".format(m3)]
    return render_template('result.html',peptide=peptide,mhc=mhc,score=score,p=p,m=m,i=i,m3=m3,path=path)




@app.route('/',methods=['GET','POST'])
def home():
    if request.method == 'POST':
        peptide = request.form['peptide']
        mhc = request.form['mhc']
        print(peptide,mhc,type(peptide),type(mhc))
        return redirect(url_for('result',peptide=peptide,mhc=mhc))
    return render_template('submit.html')





