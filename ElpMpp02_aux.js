"use strict";

function init() {
    var jd0 = 2444239.5;
    var djd = 2000.0;
    document.getElementById("jd0").value = jd0;
    document.getElementById("djd").value = djd;
    var n = parseInt(document.getElementById("n").value);
    var LLR = document.getElementById("LLR").checked;
    compute_ElpMpp02(jd0, djd, n, LLR);
}

function processForm(form) {
    var res = document.getElementById("result");
    res.innerHTML = "";
    var jd0 = parseFloat(form.jd0.value);
    var djd = parseFloat(form.djd.value);
    var n = parseInt(form.n.value);
    var LLR = form.LLR.checked;
    var message = "Invalid Julian date!";
    sanityCheck(jd0, "jd0", message, "result");
    message = "Invalid step size!";
    sanityCheck(djd, "djd", message, "result");
    if (res.innerHTML == "") {
        compute_ElpMpp02(jd0, djd, n, LLR);
    }
}

function compute_ElpMpp02(jd0, djd, n, LLR) {
    var res = document.getElementById('result');
    res.innerHTML = "Please wait...";
    var txt = "<h3>Result:</h3>"
    txt += "<table>";
    txt += "<tr><th>JD</th> <th>Date, Time<br /> (YYYY-MM-DD, HH:MM:SS)</th> <th>X (km)</th> <th>Y (km)</th> <th>Z (km)</th> </tr>";
    var T0 = (jd0 - 2451545.0)/36525.0;
    var dT = djd/36525.0;
    for (var i=0; i<n; i++) {
        var jd = jd0 + i*djd;
        var date = CalDat(jd - 2451545.0);
        var T = T0 + i*dT;
        var moon; 
        if (LLR) {
            moon = getX2000_LLR(T);
        } else {
            moon = getX2000_DE(T);
        }
        txt += "<tr> <td>"+jd.toFixed(5)+"</td><td>" + 
            date.dateString+", "+date.timeString+"</td> <td>"+
            moon.X.toFixed(5) + "</td> <td>" + 
            moon.Y.toFixed(5) + "</td> <td>" + 
            moon.Z.toFixed(5)+"</td> </tr>";
    }
    txt += "</table> <br />";
    res.innerHTML = txt;
    txt = "<p>X, Y, Z are the geocentric coordinates of the Moon with respect to the mean ecliptic and equinox of J2000.0.</p>";
    txt += "<p>Times are in barycentric dynamical time (TDB).</p>"
    res.innerHTML += txt;
}

// sanity check
// If there are errors, print message in red at the place 
// indicated by the id errid
function sanityCheck(x,inputId,message,errid) {
    var input = document.getElementById(inputId);
    input.style.backgroundColor = "white";
    if (isNaN(x)) {
        input.style.backgroundColor = "#e2a8a8";
        var text = '<p style="color:red;">'+message+'</p>';
        document.getElementById(errid).innerHTML += text;
    }
}
    
function CalDat(D) {
    var a,b,c,d,e,f;
    // Convert Julian day number to calendar date
    a = Math.floor(D + 2451545.5);
    if (a < 0) {
        return CalDatNegativeJD(D + 2451545);
    }
    if (a < 2299161) { //Julian calendar
        b = 0; c = a+1524;
    } else { // Gregorian calendar
        b = Math.floor((a-1867216.25)/36524.25);
        c = a + b - Math.floor(0.25*b) + 1525;
    }
    d = Math.floor((c-122.1)/365.25);
    if (d < 0) {d++;}
    e = 365*d + Math.floor(0.25*d);
    f = Math.floor((c-e)/30.6001);
    if (f < 0) {f++;}
    var dd = c-e - Math.floor(30.6001*f);
    var mm = f - 1 - 12*Math.floor(f/14+1e-5);
    var yy = d - 4715 - Math.floor((7+mm)/10+1e-5);
    var dateString = generateDateString(yy,mm,dd);
    var FracOfDay = D - Math.floor(D+0.5) + 0.5;
    var Hour = 24*FracOfDay;
    var h = Math.floor(Hour);
    var m = Math.floor(60*(Hour-h));
    var s = (Hour - h - m/60)*3600;
    var timeString = generateTimeString(h,m,s);
    return {yy:yy, mm:mm, dd:dd, h:h, m:m, s:s, 
           dateString:dateString, timeString:timeString};
}

//----------------------------------------------------------
// CalDat: Calendar date and time from Julian date JD with JD<0
// 
// yy,mm,dd Calendar date components
// h, m, s hour, min, sec.
// 
//-------------------------------------------------
function CalDatNegativeJD(jd) {
    var mjd = -Math.floor(jd+0.5);
    var md = mjd - Math.floor(mjd/1461);
    var dyear = Math.floor(md/(365+1e-10)) + 1;
    var yyyy = -4712 - dyear;
    var mjd0 = dyear*365 + Math.floor(dyear/4) + 1;
    var dFromY = mjd0 - mjd;
    var monthTable;
    if (dyear % 4 ==0) {
       monthTable = [0, 31, 60, 91, 121, 152, 182, 213, 244, 
                    274, 305, 335, 366];
    } else {
       monthTable = [0, 31, 59, 90, 120, 151, 181, 212, 243, 
                    273, 304, 334, 365];
    }
    var i,mm,dd;
    for (i=1; i<13; i++) {
        if (dFromY <= monthTable[i]) {
            mm = i;
            dd = dFromY - monthTable[i-1];
            break;
        }
    }
    var dateString = generateDateString(yyyy,mm,dd);
    var FracOfDay = 0.5+ (jd + mjd);
    var Hour = 24*FracOfDay;
    var h = Math.floor(Hour);
    var m = Math.floor(60*(Hour-h));
    var s = (Hour - h - m/60)*3600;
    var timeString = generateTimeString(h,m,s);
    return {yy:yyyy, mm:mm, dd:dd, h:h, m:m, s:s, 
           dateString:dateString, timeString:timeString};
}
    
// Generate date string from yyyy, mm and dd:
// return yyyy-mm-dd
function generateDateString(yyyy,mm,dd) {
    var absy = Math.abs(yyyy);
    if (absy < 10) {
        absy = "000"+absy;
    } else if (absy < 100) {
        absy = "00"+absy;
    } else if (absy < 1000) {
        absy = "0"+absy;
    } else {
        absy = absy.toString();
    }
    var yStr = absy;
    if (yyyy < 0) {yStr = "-"+yStr;}
    var mmString = mm.toString();
    if (mm < 10) {mmString = "0"+mmString;}
    var ddString = dd.toString();
    if (dd < 10) {ddString = "0"+ddString;}
    return yStr+"-"+mmString+"-"+ddString;
}

// Generate time string from h,m,s: 
// return hh:mm:ss 
function generateTimeString(h,m,s) {
    var hround = h + m/60 + (s+0.5)/3600;
    var hh = Math.floor(hround);
    var mm = Math.floor((hround-hh)*60);
    var ss= Math.floor(3600*(hround-hh-mm/60));
    hh = hh.toString(); mm = mm.toString(); ss = ss.toString();
    if (hh.length < 2) {hh = "0"+hh;}
    if (mm.length < 2) {mm = "0"+mm;}
    if (ss.length < 2) {ss = "0"+ss;}
    return hh+":"+mm+":"+ss;
}