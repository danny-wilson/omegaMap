#ifndef _MY_DOUBLE_H_
#define _MY_DOUBLE_H_

#include <limits>
#include <math.h>
#include <myerror.h>

/*	This class behaves to the user like a non-negative double, but
	is stored internally as the natural logarithm. Standard mathematical
	operations are performed on the logarithm of the number so that it
	should not underflow or overflow like a double. */
class mydouble {
protected:
	double _log;
	bool zero;
public:
	/*Default constructor*/
	mydouble() {
		zero = false;
	};
	/*Copy constructor*/
	mydouble(const double &_doub) {
		zero = false;
		if(_doub<0.0) myutils::error("mydouble::mydouble(const double&): cannot initialize with negative number");
		if(_doub==0.0) setzero();
		else _log = log(_doub);
	};
	/*Copy constructor*/
	mydouble(const mydouble &_mydoub) {
		zero = _mydoub.zero;
		_log = _mydoub._log;
	}
	/*Conversion operator*/
	operator double const() {
		return (zero) ? 0.0 : exp(_log);
	};
	/*Assignment operator*/
	mydouble& operator=(const double &_doub) {
		zero = false;
		if(_doub<0.0) myutils::error("mydouble::operator=(const double&): cannot assign a negative number");
		if(_doub==0.0) setzero();
		else _log = log(_doub);
		return *this;
	}
	/*Assignment operator*/
	mydouble& operator=(const mydouble &_mydoub) {
		zero = _mydoub.zero;
		_log = _mydoub._log;
		return *this;
	}

	void setlog(const double &log) {
		zero = false;
		_log = log;
	}
	mydouble& setzero() {
		zero = true;
		_log = -std::numeric_limits<double>::max();
		return *this;
	}
	bool iszero() {
		return zero;
	}
	
	/*** MULTIPLICATION ***/
	mydouble operator*(const double &dbl) const {
		return operator*(mydouble(dbl));
	}
	mydouble operator*(const mydouble &mydbl) const {
		mydouble a;
		if(zero || mydbl.zero) a.setzero();
		else a.setlog(_log + mydbl._log);
		return a;
	}
	mydouble& operator*=(const double &dbl) {
		if(zero || dbl==0.0) setzero();
		else _log += mydouble(dbl)._log;
		return *this;
	}
	mydouble& operator*=(const mydouble &mydbl) {
		if(zero || mydbl.zero) setzero();
		else _log += mydbl._log;
		return *this;
	}

	/*** DIVISION ***/
	mydouble operator/(const double &dbl) const {
		return operator/(mydouble(dbl));
	}
	mydouble operator/(const mydouble &mydbl) const {
		mydouble a;
		if(zero) a.setzero();
		else if(mydbl.zero) error("mydouble::operator/(const mydouble&): division by zero");
		else a.setlog(_log - mydbl._log);
		return a;
	}
	mydouble& operator/=(const double &dbl) {
		if(dbl==0.0) error("mydouble::operator/=(const double&): division by zero");
		else if(!zero) _log -= mydouble(dbl)._log;
		return *this;
	}
	mydouble& operator/=(const mydouble &mydbl) {
		if(mydbl.zero) error("mydouble::operator/=(const mydouble&): division by zero");
		else if(!zero) _log -= mydbl._log;
		return *this;
	}

	/*** ADDITION ***/
	mydouble operator+(const double &dbl) const {
		if(dbl==0.0) return mydouble(*this);
		return operator+(mydouble(dbl));
	}
	mydouble operator+(const mydouble &mydbl) const {
		mydouble a;
		if(zero) a = mydouble(mydbl);
		else if(mydbl.zero) a = mydouble(*this);
		else {
			double diff = _log - mydbl._log;
			if(diff==0.0) a.setlog(log(2.0) + _log);
			else if(diff<0.0) a.setlog(mydbl._log + log(1.0 + exp(diff)));
			else a.setlog(_log + log(1.0 + exp(-diff)));
		}
		return a;
	}
	mydouble& operator+=(const double &dbl) {
		if(dbl==0.0) return *this;
		return operator+=(mydouble(dbl));
	}
	mydouble& operator+=(const mydouble &mydbl) {
		if(zero) *this = mydbl;
		else if(!mydbl.zero) {
			double diff = _log - mydbl._log;
			if(diff==0.0) _log += log(2.0);
			else if(diff<0.0) _log = mydbl._log + log(1.0 + exp(diff));
			else _log += log(1.0 + exp(-diff));
		}
		return *this;
	}

	/*** SUBTRACTION - warning cannot have negative numbers ***/
	mydouble operator-(const double &dbl) const {
		if(dbl==0.0) return mydouble(*this);
		return operator-(mydouble(dbl));
	}
	mydouble operator-(const mydouble &mydbl) const {
		mydouble a;
		if(mydbl.zero) a = mydouble(*this);
		else if(zero) error("mydouble::operator-(const mydouble&): subtracting a positive number from zero");
		else {
			/* diff must always be positive */
			double diff = _log - mydbl._log;
			if(diff==0.0) a.setzero();
			else if(diff<0.0) myutils::error("mydouble::operator-(const mydouble&) cannot handle negative numbers");
			else a.setlog(_log + log(1.0 - exp(-diff)));
		}
		return a;
	}
	mydouble& operator-=(const double &dbl) {
		if(dbl==0.0) return *this;
		return operator-=(mydouble(dbl));
	}
	mydouble& operator-=(const mydouble &mydbl) {
		if(!mydbl.zero) {
			if(zero) error("mydouble::operator-=(const mydouble&): subtracting a positive number from zero");
			/* diff must always be positive */
			double diff = _log - mydbl._log;
			if(diff==0.0) setzero();
			else if(diff<0.0) myutils::error("mydouble::operator-=(const mydouble&) cannot handle negative numbers");
			else _log += log(1.0 - exp(-diff));
		}
		return *this;
	}

	/*** SPECIAL OPERATIONS ***/
	double LOG() const {
		return _log;
	}
	/* Caution: ^ has lower precedence than /*+- */
	mydouble operator^(const double &dbl) const {
		mydouble a;
		if(zero) a.setzero();
		else a.setlog(_log * dbl);
		return a;
	}
	/* Caution: ^ has lower precedence than /*+- */
	mydouble operator^(const mydouble &mydbl) const {
		mydouble a;
		if(zero) a.setzero();
		else a.setlog(_log * exp(mydbl._log));
		return a;
	}
	mydouble& operator^=(const double &dbl) {
		if(!zero) _log *= dbl;
		return *this;
	}
	mydouble& operator^=(const mydouble &mydbl) {
		if(!zero) _log *= exp(mydbl._log);
		return *this;
	}
};

#endif//_MY_DOUBLE_H_
