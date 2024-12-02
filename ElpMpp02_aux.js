"use strict";

function init(v=0) {
    let jd0 = 2444239.5;
    let djd = 2000.0;
    document.getElementById("jd0").value = jd0;
    document.getElementById("djd").value = djd;
    let n = parseInt(document.getElementById("n").value, 10);
    let LLR = document.getElementById("LLR").checked;
    compute_ElpMpp02(jd0, djd, n, LLR, v);
}

function processForm(form, v=0) {
    let res = document.getElementById("result");
    res.innerHTML = "";
    let jd0 = parseFloat(form.jd0.value);
    let djd = parseFloat(form.djd.value);
    let n = parseInt(form.n.value);
    let LLR = form.LLR.checked;
    let message = "Invalid Julian date!";
    sanityCheck(jd0, "jd0", message, "result");
    message = "Invalid step size!";
    sanityCheck(djd, "djd", message, "result");
    if (res.innerHTML == "") {
        compute_ElpMpp02(jd0, djd, n, LLR, v);
    }
}

function compute_ElpMpp02(jd0, djd, n, LLR, v) {
    let res = document.getElementById('result');
    res.innerHTML = "Please wait...";
    let txt = "<h3>Result:</h3>"
    txt += "<table>";
    if (v==0) {
        txt += "<tr><th>JD</th> <th>Date, Time<br /> (YYYY-MM-DD, HH:MM:SS)</th> <th>X (km)</th> <th>Y (km)</th> <th>Z (km)</th> </tr>";
    } else {
        txt += "<tr><th>JD</th> <th>Date, Time<br /> (YYYY-MM-DD, HH:MM:SS)</th> <th>X (km)<br />Vx (km/day)</th> <th>Y (km)<br />Vy (km/day)</th> <th>Z (km)<br />Vz (km/day)</th> </tr>";
    }
    let T0 = (jd0 - 2451545.0)/36525.0;
    let dT = djd/36525.0;
    for (let i=0; i<n; i++) {
        let jd = jd0 + i*djd;
        let date = CalDat(jd - 2451545.0);
        let T = T0 + i*dT;
        let moon; 
        if (LLR) {
            moon = (v==0 ? getX2000_LLR(T):getX2000_Xdot2000_LLRm(T));
        } else {
            moon = (v==0 ? getX2000_DE(T):getX2000_Xdot2000_DEm(T));
        }
        txt += '<tr> <td>'+jd.toFixed(5)+'</td><td>' + 
            date.dateString+', '+date.timeString+'</td> <td>'+
            moon.X.toFixed(5) + (v==0 ? '</td>':'<br />'+moon.Xdot.toFixed(5)+'</td>') + '<td>' + 
            moon.Y.toFixed(5) + (v==0 ? '</td>':'<br />'+moon.Ydot.toFixed(5)+'</td>') + '<td>' + 
            moon.Z.toFixed(5) + (v==0 ? '</td>':'<br />'+moon.Zdot.toFixed(5)+'</td>') + '</tr>';
    }
    txt += '</table> <br />';
    txt += '<p>X, Y, Z are the geocentric coordinates of the Moon with respect to the mean ecliptic and equinox of J2000.0'; 
    txt += (v==0 ? '.':"; Vx, Vy, Vz are components of Moon's velocity.") + '</p>';
    txt += '<p>Times are in barycentric dynamical time (TDB).</p>';
    res.innerHTML = txt;
}

// sanity check
// If there are errors, print message in red at the place 
// indicated by the id errid
function sanityCheck(x,inputId,message,errid) {
    let input = document.getElementById(inputId);
    input.style.backgroundColor = "transparent";
    if (isNaN(x)) {
        input.style.backgroundColor = "#e2a8a8";
        let text = '<p style="color:red;">'+message+'</p>';
        document.getElementById(errid).innerHTML += text;
    }
}
    
function CalDat(D) {
    let a,b,c,d,e,f;
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
    let dd = c-e - Math.floor(30.6001*f);
    let mm = f - 1 - 12*Math.floor(f/14+1e-5);
    let yy = d - 4715 - Math.floor((7+mm)/10+1e-5);
    let dateString = generateDateString(yy,mm,dd);
    let FracOfDay = D - Math.floor(D+0.5) + 0.5;
    let Hour = 24*FracOfDay;
    let h = Math.floor(Hour);
    let m = Math.floor(60*(Hour-h));
    let s = (Hour - h - m/60)*3600;
    let timeString = generateTimeString(h,m,s);
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
    let mjd = -Math.floor(jd+0.5);
    let md = mjd - Math.floor(mjd/1461);
    let dyear = Math.floor(md/(365+1e-10)) + 1;
    let yyyy = -4712 - dyear;
    let mjd0 = dyear*365 + Math.floor(dyear/4) + 1;
    let dFromY = mjd0 - mjd;
    let monthTable;
    if (dyear % 4 ==0) {
       monthTable = [0, 31, 60, 91, 121, 152, 182, 213, 244, 
                    274, 305, 335, 366];
    } else {
       monthTable = [0, 31, 59, 90, 120, 151, 181, 212, 243, 
                    273, 304, 334, 365];
    }
    let i,mm,dd;
    for (i=1; i<13; i++) {
        if (dFromY <= monthTable[i]) {
            mm = i;
            dd = dFromY - monthTable[i-1];
            break;
        }
    }
    let dateString = generateDateString(yyyy,mm,dd);
    let FracOfDay = 0.5+ (jd + mjd);
    let Hour = 24*FracOfDay;
    let h = Math.floor(Hour);
    let m = Math.floor(60*(Hour-h));
    let s = (Hour - h - m/60)*3600;
    let timeString = generateTimeString(h,m,s);
    return {yy:yyyy, mm:mm, dd:dd, h:h, m:m, s:s, 
           dateString:dateString, timeString:timeString};
}
    
// Generate date string from yyyy, mm and dd:
// return yyyy-mm-dd
function generateDateString(yyyy,mm,dd) {
    let absy = Math.abs(yyyy);
    if (absy < 10) {
        absy = "000"+absy;
    } else if (absy < 100) {
        absy = "00"+absy;
    } else if (absy < 1000) {
        absy = "0"+absy;
    } else {
        absy = absy.toString();
    }
    let yStr = absy;
    if (yyyy < 0) {yStr = "-"+yStr;}
    let mmString = mm.toString();
    if (mm < 10) {mmString = "0"+mmString;}
    let ddString = dd.toString();
    if (dd < 10) {ddString = "0"+ddString;}
    return yStr+"-"+mmString+"-"+ddString;
}

// Generate time string from h,m,s: 
// return hh:mm:ss 
function generateTimeString(h,m,s) {
    let hround = h + m/60 + (s+0.5)/3600;
    let hh = Math.floor(hround);
    let mm = Math.floor((hround-hh)*60);
    let ss= Math.floor(3600*(hround-hh-mm/60));
    hh = hh.toString(); mm = mm.toString(); ss = ss.toString();
    if (hh.length < 2) {hh = "0"+hh;}
    if (mm.length < 2) {mm = "0"+mm;}
    if (ss.length < 2) {ss = "0"+ss;}
    return hh+":"+mm+":"+ss;
}